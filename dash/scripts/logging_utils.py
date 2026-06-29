"""
Shared logging helpers for standalone data preparation scripts.

Each script configures the root logger once and proxies legacy `print`
calls so that every console message also follows the unified logging
format (timestamp + severity). The proxy keeps the original call
signature (`print("msg", sep=" ", end="\\n")`) so that existing scripts
do not need invasive rewrites.
"""

from __future__ import annotations

import logging
from typing import Callable

DEFAULT_FORMAT = "%(asctime)s [%(levelname)s] %(message)s"


def configure_script_logger(name: str, level: int = logging.INFO) -> logging.Logger:
    """
    Configure logging for CLI-style scripts.

    The first call sets up the root logger with a console handler; subsequent
    calls simply return module-specific loggers to avoid duplicate handlers.
    """
    root_logger = logging.getLogger()
    if not root_logger.handlers:
        logging.basicConfig(level=level, format=DEFAULT_FORMAT)
    logger = logging.getLogger(name)
    logger.setLevel(level)
    return logger


def create_print_proxy(logger: logging.Logger) -> Callable[..., None]:
    """
    Return a callable that mimics `print()` but writes through logging.

    Scripts can simply assign `print = create_print_proxy(logger)` at import
    time, and every old `print()` call will show up in log files/STDOUT with
    consistent formatting.
    """

    def _proxy(*args, level: int | None = None, **kwargs) -> None:
        sep = kwargs.pop("sep", " ")
        end = kwargs.pop("end", "\n")
        kwargs.pop("file", None)
        kwargs.pop("flush", None)

        message = sep.join(str(arg) for arg in args)
        if end and end != "\n":
            message = f"{message}{end}"

        normalized = message.lstrip()
        derived_level = level or _guess_level(normalized)
        logger.log(derived_level, message)

    return _proxy


ERROR_PREFIXES = ("\N{CROSS MARK}", "error", "failed", "traceback")
WARNING_PREFIXES = ("\N{WARNING SIGN}", "warning", "warn")
INFO_PREFIXES = ("\N{WHITE HEAVY CHECK MARK}", "\N{PARTY POPPER}")


def _guess_level(message: str) -> int:
    """Map emoji/prefixes to a reasonable logging level."""
    lowered = message.lower()
    if lowered.startswith(ERROR_PREFIXES):
        return logging.ERROR
    if lowered.startswith(WARNING_PREFIXES):
        return logging.WARNING
    if lowered.startswith(INFO_PREFIXES):
        return logging.INFO
    return logging.INFO
