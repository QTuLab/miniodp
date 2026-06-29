import logging
import os
import signal
import time
from pathlib import Path
from typing import Optional

import psutil


logger = logging.getLogger("miniodp.memory_watchdog")
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s [memory_watchdog] %(message)s",
)


def read_int(path: Path) -> Optional[int]:
    try:
        return int(path.read_text().strip())
    except FileNotFoundError:
        return None
    except Exception as exc:  # noqa: BLE001
        logger.warning("Failed to read %s: %s", path, exc)
        return None


def read_cgroup_usage_bytes() -> Optional[int]:
    """Read container memory usage from cgroup when available."""
    for path in (
        Path("/sys/fs/cgroup/memory.current"),  # cgroup v2
        Path("/sys/fs/cgroup/memory/memory.usage_in_bytes"),  # cgroup v1
    ):
        value = read_int(path)
        if value is not None:
            return value
    return None


def read_process_tree_usage_bytes(master_pid: int) -> Optional[int]:
    try:
        master = psutil.Process(master_pid)
    except psutil.NoSuchProcess:
        return None
    try:
        total = master.memory_info().rss
        for child in master.children(recursive=True):
            try:
                total += child.memory_info().rss
            except psutil.NoSuchProcess:
                continue
        return total
    except psutil.NoSuchProcess:
        return None


def get_used_bytes(master_pid: int) -> Optional[int]:
    usage = read_cgroup_usage_bytes()
    if usage is not None:
        return usage

    usage = read_process_tree_usage_bytes(master_pid)
    if usage is not None:
        return usage

    try:
        return psutil.virtual_memory().used
    except Exception as exc:  # noqa: BLE001
        logger.warning("Failed to read memory usage: %s", exc)
        return None


def load_master_pid(pid_file: Path) -> Optional[int]:
    pid = read_int(pid_file)
    if pid is None:
        return None
    try:
        proc = psutil.Process(pid)
        if proc.is_running():
            return pid
    except psutil.NoSuchProcess:
        return None
    return None


def to_bytes(gb: float) -> int:
    return int(gb * 1024 * 1024 * 1024)


def main() -> None:
    threshold_gb = float(os.environ.get("MEM_WATCHDOG_THRESHOLD_GB", "40"))
    interval_sec = int(os.environ.get("MEM_WATCHDOG_INTERVAL_SEC", "15"))
    trigger_count_needed = int(os.environ.get("MEM_WATCHDOG_TRIGGER_COUNT", "3"))
    cooldown_sec = int(os.environ.get("MEM_WATCHDOG_COOLDOWN_SEC", "300"))
    pid_file = Path(os.environ.get("GUNICORN_PID_FILE", "/app/logs/gunicorn.pid"))

    threshold_bytes = to_bytes(threshold_gb)

    logger.info(
        "Starting memory watchdog: threshold %.1f GB, interval %ss, trigger count %s, cooldown %ss, pid file %s",
        threshold_gb,
        interval_sec,
        trigger_count_needed,
        cooldown_sec,
        pid_file,
    )

    strikes = 0
    cooldown_until = 0.0

    while True:
        now = time.time()
        if now < cooldown_until:
            time.sleep(interval_sec)
            continue

        master_pid = load_master_pid(pid_file)
        if master_pid is None:
            logger.info("Gunicorn master pid not found yet; retrying later")
            time.sleep(interval_sec)
            continue

        used = get_used_bytes(master_pid)
        if used is None:
            logger.info("Memory usage is unavailable; retrying later")
            time.sleep(interval_sec)
            continue

        used_gb = used / (1024 * 1024 * 1024)

        if used >= threshold_bytes:
            strikes += 1
            logger.warning(
                "Memory threshold exceeded: current %.2f GB, strike %d/%d",
                used_gb,
                strikes,
                trigger_count_needed,
            )
            if strikes >= trigger_count_needed:
                logger.warning(
                    "Trigger reached; sending HUP to gunicorn master (%d) for graceful worker restart",
                    master_pid,
                )
                try:
                    os.kill(master_pid, signal.SIGHUP)
                except ProcessLookupError:
                    logger.warning("Failed to send signal; gunicorn master may have exited")
                strikes = 0
                cooldown_until = time.time() + cooldown_sec
        else:
            strikes = 0

        time.sleep(interval_sec)


if __name__ == "__main__":
    main()
