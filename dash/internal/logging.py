import logging
import os
from pathlib import Path


def configure_logging() -> logging.Logger:
    """Configure logging outputs and return module logger."""
    log_level_name = os.environ.get("DASH_LOG_LEVEL", "INFO").upper()
    log_format = "%(asctime)s [%(levelname)s] %(name)s: %(message)s"
    log_level = getattr(logging, log_level_name, logging.INFO)
    root_logger = logging.getLogger()

    # In production (e.g., Gunicorn), the process typically pre-configures root handlers.
    # Avoid overriding them here to prevent duplicated/misrouted logs.
    if not root_logger.handlers:
        handlers = [logging.StreamHandler()]
        log_file = os.environ.get("DASH_LOG_FILE")
        if log_file:
            log_path = Path(log_file)
            log_path.parent.mkdir(parents=True, exist_ok=True)
            handlers.append(logging.FileHandler(log_path))
        logging.basicConfig(level=log_level, format=log_format, handlers=handlers)
    else:
        root_logger.setLevel(log_level)

    return logging.getLogger("miniodp.dash")

