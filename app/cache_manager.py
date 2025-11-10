# cache_manager.py
from __future__ import annotations
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

class CacheGate:
    """
    Pause prefetch when cache usage >= high; resume when usage <= low.
    Pure polling. No eviction; processing must delete .sra files.
    """
    def __init__(self, cache_dir: Path, high_bytes: int, low_bytes: int, poll_secs: float = 3.0):
        assert low_bytes < high_bytes, "low watermark must be < high watermark"
        self.cache_dir = cache_dir.expanduser().resolve()
        self.high = int(high_bytes)
        self.low = int(low_bytes)
        self.poll = float(poll_secs)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self._fh = None
        self._lockfile = self.cache_dir / LOCKFILE

    def _lock(self):
        if self._fh is None:
            self._fh = open(self._lockfile, "a+")
        fcntl.flock(self._fh.fileno(), fcntl.LOCK_EX)

    def _unlock(self):
        if self._fh:
            fcntl.flock(self._fh.fileno(), fcntl.LOCK_UN)

    def wait_for_window(self, log_fn, tag: str):
        """Block if usage >= high; resume when usage <= low."""
        gb = 1024**3
        while True:
            self._lock()
            try:
                used = cache_usage_bytes(self.cache_dir)
                if used < self.high:
                    return
                log_fn(f"[{tag}] Prefetch paused: cache {used/gb:.1f} GiB ≥ high {self.high/gb:.1f} GiB. "
                       f"Waiting until ≤ low {self.low/gb:.1f} GiB …")
                self._unlock()
                while True:
                    time.sleep(self.poll)
                    used = cache_usage_bytes(self.cache_dir)
                    if used <= self.low:
                        return
            finally:
                self._unlock()
