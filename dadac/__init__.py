import platform
if platform.platform() != "Linux":
    raise(ModuleNotFoundError(f"Implementation actually supports only Linux, you are running on {platform.platform()} "))
from .dadac import *
