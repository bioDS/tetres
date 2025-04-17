"""
Module contains decorators used across the package.
"""
from typing import get_args, Literal
from functools import wraps


def validate_literal_args(**literal_checks):
    """
    Decorator to validate function arguments against Literal options.

    Example:
        @validate_literal_args(mds_type=_MDS_TYPES, dist_type=_DIST)
    """
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            for arg_name, literal_type in literal_checks.items():
                value = kwargs.get(arg_name)
                if value is not None:
                    options = get_args(literal_type)
                    if value not in options:
                        raise ValueError(
                            f"Invalid value for '{arg_name}': '{value}'. "
                            f"Choose from: {', '.join(map(str, options))}"
                        )
            return func(*args, **kwargs)
        return wrapper
    return decorator
