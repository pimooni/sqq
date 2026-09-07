"""Structure writers and output-directory lifecycle helpers.

Concrete services are imported from their modules so rendering and cleanup do
not acquire a circular package-initialization dependency.
"""

__all__: list[str] = []
