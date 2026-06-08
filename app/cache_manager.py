# cache_manager.py
from __future__ import annotations
import os
import time
import fcntl
from pathlib import Path

LOCKFILE = ".hulk_cache.lock"


def cache_usage_bytes(cache_dir: Path) -> int:
    """
    Calculates the total storage space currently occupied by SRA files within the designated cache directory.

    This function recursively aggregates the size of all `.sra` files directly in the root
    or located within standard SRR subdirectories.

    Args:
        cache_dir (Path): The root directory of the SRA cache.

    Returns:
        int: The total size of the cached files in bytes.
    """
    cache_dir = cache_dir.expanduser().resolve()
    total = 0

    # Process files directly in the cache root
    for p in cache_dir.glob("*.sra"):
        try:
            total += p.stat().st_size
        except FileNotFoundError:
            pass

    # Process files within standard 'SRR*' subdirectories
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
    """
    Determines the available free space on the filesystem hosting the provided path.

    Args:
        path (Path): A path located on the target filesystem.

    Returns:
        int | None: The available free space in bytes, or None if the operation fails (e.g., due to OS restrictions).
    """
    try:
        path = path.expanduser().resolve()
        st = os.statvfs(str(path))
        return st.f_bavail * st.f_frsize
    except OSError:
        return None


class CacheGate:
    """
    Manages cache capacity by enforcing high and low watermarks to pause or resume prefetching.

    The CacheGate logic ensures that the system does not exhaust disk space. If the storage usage
    meets or exceeds the high watermark, prefetching operations will be blocked until the usage
    drops back down to the low watermark.

    Constraints managed internally:
      - Requires at least 20 GiB of free space to remain enabled.
      - Automatically caps the high watermark at (free space - 20 GiB buffer).
      - Forces the low watermark to strictly remain below the high watermark (defaults to 80% if invalid).
    """

    def __init__(self, cache_dir: Path, high_bytes: int, low_bytes: int, poll_secs: float = 3.0):
        """
        Initializes the CacheGate and validates the environmental capacity constraints.

        Args:
            cache_dir (Path): The path to the cache directory.
            high_bytes (int): The upper limit of allowed cache usage in bytes.
            low_bytes (int): The lower threshold to resume processing in bytes.
            poll_secs (float, optional): Interval in seconds to check usage when blocked. Defaults to 3.0.
        """
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

        # Disable caching if thresholds are logically invalid.
        if hb <= 0 or lb <= 0:
            self.disabled = True
            self._init_msg = "[cache] disabled: non-positive watermarks"
            return

        if free is not None:
            # Immediate shutdown if the system lacks the minimum buffer.
            if free <= buf:
                self.disabled = True
                self._init_msg = (
                    f"[cache] disabled: only {free / gb:.1f} GiB free "
                    f"(buffer={buf / gb:.0f} GiB)"
                )
                return

            # Apply capacity ceiling
            if hb >= free:
                effective_high = max(0, free - buf)
                if effective_high <= 0:
                    self.disabled = True
                    self._init_msg = (
                        f"[cache] disabled: requested high too large and "
                        f"free space {free / gb:.1f} GiB ≈ buffer"
                    )
                    return
                hb = effective_high

            # Correct logical inversions between high and low thresholds.
            if lb >= hb:
                lb = max(1, int(hb * 0.8))

        if lb >= hb:
            self.disabled = True
            self._init_msg = "[cache] disabled: low watermark >= high watermark after adjustment"
            return

        self.high = hb
        self.low = lb

    def _lock(self):
        """Acquires an exclusive file lock to ensure atomic read operations across processes."""
        if self._fh is None:
            self._fh = open(self._lockfile, "a+")
        fcntl.flock(self._fh.fileno(), fcntl.LOCK_EX)

    def _unlock(self):
        """Releases the exclusive file lock."""
        if self._fh:
            fcntl.flock(self._fh.fileno(), fcntl.LOCK_UN)

    def wait_for_window(self, log_fn, tag: str):
        """
        Blocks the calling thread if the cache size has exceeded the high watermark, releasing
        only when it falls below the low watermark.

        If the gate was initialized as disabled (due to low space or misconfiguration),
        this function operates as a non-blocking pass-through.

        Args:
            log_fn (callable): A logging function to output messages.
            tag (str): Contextual identifier used for log formatting.
        """
        gb = 1024 ** 3

        if self.disabled:
            if self._init_msg:
                log_fn(f"[{tag}] {self._init_msg}")
                self._init_msg = None
            return

        # Polling loop to block prefetching threads
        while True:
            self._lock()
            try:
                used = cache_usage_bytes(self.cache_dir)
                if used < self.high:
                    return
                log_fn(
                    f"[{tag}] Prefetch paused: cache {used / gb:.1f} GiB ≥ high {self.high / gb:.1f} GiB. "
                    f"Waiting until ≤ low {self.low / gb:.1f} GiB …"
                )
                self._unlock()

                # Sleep and verify loop while cache clears
                while True:
                    time.sleep(self.poll)
                    used = cache_usage_bytes(self.cache_dir)
                    if used <= self.low:
                        return
            finally:
                self._unlock()