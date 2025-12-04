# cache_manager.py
from __future__ import annotations
import os
import time
import fcntl
from pathlib import Path

LOCKFILE = ".hulk_cache.lock"

def cache_usage_bytes(cache_dir: Path) -> int:
    """Sum sizes of *.sra in cache_dir (handles SRR/SRR.sra subdirs too)."""
    cache_dir = cache_dir.expanduser().resolve()
    total = 0
    for p in cache_dir.glob("*.sra"):
        try:
            total += p.stat().st_size
        except FileNotFoundError:
            pass
    for d in cache_dir.glob("SRR*"):
        if d.is_dir():
            sra = d / f"{d.name}.sra"
            try:
                if sra.exists():
                    total += sra.stat().st_size
            except FileNotFoundError:
                pass
    return total

def _free_space_bytes(path: Path) -> int | None:
    """Return free space in bytes for the filesystem containing `path`."""
    try:
        path = path.expanduser().resolve()
        st = os.statvfs(str(path))
        return st.f_bavail * st.f_frsize
    except OSError:
        return None

class CacheGate:
    """
    Pause prefetch when cache usage >= high; resume when usage <= low.

    All *available space* and watermark sanity is handled here:
      • If free space <= 20 GiB → gate is disabled (no caching).
      • If requested high >= free space → cap high at (free - 20 GiB).
      • If low >= high → set low to 80% of high.
      • If resulting high/low are invalid → gate is disabled.
    """
    def __init__(self, cache_dir: Path, high_bytes: int, low_bytes: int, poll_secs: float = 3.0):
        self.cache_dir = cache_dir.expanduser().resolve()
        self.poll = float(poll_secs)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self._fh = None
        self._lockfile = self.cache_dir / LOCKFILE
        self.disabled = False
        self._init_msg: str | None = None

        gb = 1024 ** 3
        free = _free_space_bytes(self.cache_dir)
        buf = 20 * gb

        hb = int(high_bytes)
        lb = int(low_bytes)

        # If requested watermarks are non-positive, just disable the gate
        if hb <= 0 or lb <= 0:
            self.disabled = True
            self._init_msg = "[cache] disabled: non-positive watermarks"
            return

        if free is not None:
            if free <= buf:
                # Not enough free space to safely cache
                self.disabled = True
                self._init_msg = (
                    f"[cache] disabled: only {free/gb:.1f} GiB free "
                    f"(buffer={buf/gb:.0f} GiB)"
                )
                return

            # Cap high at (free - buffer)
            if hb >= free:
                effective_high = max(0, free - buf)
                if effective_high <= 0:
                    self.disabled = True
                    self._init_msg = (
                        f"[cache] disabled: requested high too large and "
                        f"free space {free/gb:.1f} GiB ≈ buffer"
                    )
                    return
                hb = effective_high

            # Ensure low < high; default low to 80% of high if invalid
            if lb >= hb:
                lb = max(1, int(hb * 0.8))

        # Final sanity
        if lb >= hb:
            self.disabled = True
            self._init_msg = "[cache] disabled: low watermark >= high watermark after adjustment"
            return

        self.high = hb
        self.low = lb

    def _lock(self):
        if self._fh is None:
            self._fh = open(self._lockfile, "a+")
        fcntl.flock(self._fh.fileno(), fcntl.LOCK_EX)

    def _unlock(self):
        if self._fh:
            fcntl.flock(self._fh.fileno(), fcntl.LOCK_UN)

    def wait_for_window(self, log_fn, tag: str):
        """
        Block if usage >= high; resume when usage <= low.

        If the gate is disabled (no-cache or not enough space), this is a no-op,
        but we log the reason once via `log_fn`.
        """
        gb = 1024 ** 3

        if self.disabled:
            if self._init_msg:
                log_fn(f"[{tag}] {self._init_msg}")
                self._init_msg = None
            return

        while True:
            self._lock()
            try:
                used = cache_usage_bytes(self.cache_dir)
                if used < self.high:
                    return
                log_fn(
                    f"[{tag}] Prefetch paused: cache {used/gb:.1f} GiB ≥ high {self.high/gb:.1f} GiB. "
                    f"Waiting until ≤ low {self.low/gb:.1f} GiB …"
                )
                self._unlock()
                while True:
                    time.sleep(self.poll)
                    used = cache_usage_bytes(self.cache_dir)
                    if used <= self.low:
                        return
            finally:
                self._unlock()
