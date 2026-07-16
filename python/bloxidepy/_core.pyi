from pathlib import Path

class Config:
    R: float
    gamma: float
    Pr: float
    p_e: float
    u_e: float
    T_e: float
    T_wall: float
    x: float

    def __init__(
        self,
        R: float,
        gamma: float,
        Pr: float,
        p_e: float,
        u_e: float,
        T_e: float,
        T_wall: float,
        x: float,
    ) -> None: ...
    def __eq__(self, other: object) -> bool: ...
    def __repr__(self) -> str: ...

def get_heat_transfer(x: float) -> float: ...
def read_config_file(path: str | Path) -> Config: ...
