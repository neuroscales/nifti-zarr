"""Converts units between zarr and other specifications."""

ome_valid_units = {
    "space": [
        "angstrom",
        "attometer",
        "centimeter",
        "decimeter",
        "exameter",
        "femtometer",
        "foot",
        "gigameter",
        "hectometer",
        "inch",
        "kilometer",
        "megameter",
        "meter",
        "micrometer",
        "mile",
        "millimeter",
        "nanometer",
        "parsec",
        "petameter",
        "picometer",
        "terameter",
        "yard",
        "yoctometer",
        "yottameter",
        "zeptometer",
        "zettameter",
    ],
    "time": [
        "attosecond",
        "centisecond",
        "day",
        "decisecond",
        "exasecond",
        "femtosecond",
        "gigasecond",
        "hectosecond",
        "hour",
        "kilosecond",
        "megasecond",
        "microsecond",
        "millisecond",
        "minute",
        "nanosecond",
        "petasecond",
        "picosecond",
        "second",
        "terasecond",
        "yoctosecond",
        "yottasecond",
        "zeptosecond",
        "zettasecond",
    ],
}

nifti_valid_units = [
    "unknown",
    "meter",
    "mm",
    "micron",
    "sec",
    "msec",
    "usec",
    "hz",
    "ppm",
    "rads",
]

si_prefix_short2long = {
    "Q": "quetta",
    "R": "ronna",
    "Y": "yotta",
    "Z": "zetta",
    "E": "exa",
    "P": "peta",
    "T": "tera",
    "G": "giga",
    "M": "mega",
    "K": "kilo",
    "k": "kilo",
    "H": "hecto",
    "h": "hecto",
    "D": "deca",
    "da": "deca",
    "d": "deci",
    "c": "centi",
    "m": "milli",
    "u": "micro",
    "μ": "micro",
    "n": "nano",
    "p": "pico",
    "f": "femto",
    "a": "atto",
    "z": "zepto",
    "y": "yocto",
    "r": "ronto",
    "q": "quecto",
}

si_prefix_long2short = {long: short for short, long in si_prefix_short2long.items()}

si_prefix_exponent = {
    "Q": 30,
    "R": 27,
    "Y": 24,
    "Z": 21,
    "E": 18,
    "P": 15,
    "T": 12,
    "G": 9,
    "M": 6,
    "K": 3,
    "k": 3,
    "H": 2,
    "h": 2,
    "D": 1,
    "da": 1,
    "": 0,
    "d": -1,
    "c": -2,
    "m": -3,
    "u": -6,
    "μ": -6,
    "n": -9,
    "p": -12,
    "f": -15,
    "a": -18,
    "z": -21,
    "y": -24,
    "r": -27,
    "q": -30,
}

unit_space_short2long = {
    short + "m": long + "meter" for short, long in si_prefix_short2long.items()
}
unit_space_short2long.update(
    {
        "m": "meter",
        "mi": "mile",
        "yd": "yard",
        "ft": "foot",
        "in": "inch",
        "'": "foot",
        '"': "inch",
        "Å": "angstrom",
        "pc": "parsec",
    }
)
unit_space_long2short = {long: short for short, long in unit_space_short2long.items()}
unit_space_long2short["micron"] = "u"

unit_time_short2long = {
    short + "s": long + "second" for short, long in si_prefix_short2long.items()
}
unit_time_short2long.update(
    {
        "y": "year",
        "d": "day",
        "h": "hour",
        "m": "minute",
        "s": "second",
    }
)
unit_time_long2short = {long: short for short, long in unit_time_short2long.items()}

unit_space_scale = {
    prefix + "m": 10 ** exponent for prefix, exponent in si_prefix_exponent.items()
}
unit_space_scale.update(
    {
        "mi": 1609.344,
        "yd": 0.9144,
        "ft": 0.3048,
        "'": 0.3048,
        "in": 25.4e-3,
        '"': 25.4e-3,
        "Å": 1e-10,
        "pc": 3.0857e16,
    }
)

unit_time_scale = {
    prefix + "s": 10 ** exponent for prefix, exponent in si_prefix_exponent.items()
}
unit_time_scale.update(
    {
        "y": 365.25 * 24 * 60 * 60,
        "d": 24 * 60 * 60,
        "h": 60 * 60,
        "m": 60,
    }
)


def convert_unit(value: float, src: str, dst: str) -> float:
    """Convert unit for a value."""
    src = unit_to_scale(src)
    dst = unit_to_scale(dst)
    return value * (src / dst)


def to_ome_unit(unit: str) -> str:
    """Convert unit to ome-zarr spec."""
    if unit in unit_space_short2long:
        unit = unit_space_short2long[unit]
    elif unit in unit_time_short2long:
        unit = unit_time_short2long[unit]
    elif unit in si_prefix_short2long:
        unit = si_prefix_short2long[unit]
    if unit not in (*ome_valid_units["space"], *ome_valid_units["time"]):
        raise ValueError("Unknow unit")
    return unit


def to_nifti_unit(unit: str) -> str:
    """Convert unit to nifti spec."""
    unit = to_ome_unit(unit)
    return {
        "meter": "meter",
        "millimeter": "mm",
        "micrometer": "micron",
        "second": "sec",
        "millisecond": "msec",
        "microsecond": "usec",
    }.get(unit, "unknown")


def unit_to_scale(unit: str) -> float:
    """Convert unit to scale."""
    if unit in unit_space_long2short:
        unit = unit_space_long2short[unit]
    elif unit in unit_time_long2short:
        unit = unit_time_long2short[unit]
    elif unit in si_prefix_long2short:
        unit = si_prefix_long2short[unit]
    if unit in unit_space_scale:
        unit = unit_space_scale[unit]
    elif unit in unit_time_scale:
        unit = unit_time_scale[unit]
    elif unit in si_prefix_exponent:
        unit = 10 ** si_prefix_exponent[unit]
    if isinstance(unit, str):
        raise ValueError("Unknown unit", unit)
    return unit
