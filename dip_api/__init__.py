from .penetration import PenetrationAnalisys
from .run import dip_frame
import os
import shutil

if os.name != "nt":
    target = f"{os.path.dirname(os.path.realpath(__file__))}"+"/data/dip.exe"
    os.system(f"chmod +x {target}")
    shutil.copy(target, os.getcwd())

__all__= ["dip_frame", "PenetrationAnalisys"]
