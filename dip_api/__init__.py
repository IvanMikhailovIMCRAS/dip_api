from .penetration import PenetrationAnalisys
from .run import dip_frame
import os

if os.name != "nt":
    target = f"{os.path.dirname(os.path.realpath(__file__))}"[:-7]+"/dip.exe"
    os.system(f"chmod +x {target}")

__all__= ["dip_frame", "PenetrationAnalisys"]
