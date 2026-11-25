"""Operating calendar: business days, holidays and day/night periods.

The calendar is the deterministic backbone of the schedule generators.
Times are expressed in days from 00:00 local time on 1 January of the
scenario year; integer parts index calendar days.
"""

from __future__ import annotations

import datetime as dt
from dataclasses import dataclass, field

import holidays
import numpy as np
from numpy.typing import NDArray

FloatArray = NDArray[np.float64]
IntArray = NDArray[np.int64]
BoolArray = NDArray[np.bool_]

# Day/night boundaries used by the FAA notice of intent (91 FR 65290):
# day 07:00-22:00, night 22:00-07:00.
NIGHT_START_H = 22.0
NIGHT_END_H = 7.0
DAY_HOURS = NIGHT_START_H - NIGHT_END_H  # 15 h
NIGHT_HOURS = 24.0 - DAY_HOURS  # 9 h

# Holidays on which the 2022 PEA mitigation measures preclude State Highway 4
# access restrictions (2022 PEA section 2.1.3.5). Dates are resolved per year.
PEA_RESTRICTED_HOLIDAYS = (
    "Memorial Day",
    "Labor Day",
    "Independence Day",
    "Martin Luther King Jr. Day",
    "Washington's Birthday",
    "Texas Independence Day",
    "Cesar Chavez Day",
    "Emancipation Day In Texas",
    "Juneteenth National Independence Day",
)


@dataclass(frozen=True)
class OpCalendar:
    """Calendar for one scenario year.

    Attributes
    ----------
    year:
        Calendar year used to resolve weekdays and holidays.
    n_days:
        365 or 366.
    business:
        Boolean mask over calendar days: Monday-Friday and not a U.S.
        federal holiday.
    weekend:
        Boolean mask for Saturday and Sunday.
    holiday:
        Boolean mask for U.S. federal holidays (observed).
    summer:
        Boolean mask for days from Memorial Day through Labor Day inclusive,
        the period in which the 2022 PEA mitigation measures preclude
        weekend access restrictions.
    """

    year: int
    n_days: int
    business: BoolArray = field(repr=False)
    weekend: BoolArray = field(repr=False)
    holiday: BoolArray = field(repr=False)
    summer: BoolArray = field(repr=False)
    protected_holiday: BoolArray = field(repr=False)

    @property
    def business_days(self) -> IntArray:
        """Indices of business days."""
        return np.flatnonzero(self.business).astype(np.int64)

    @property
    def n_business(self) -> int:
        return int(self.business.sum())

    def day_of_week(self) -> IntArray:
        start = dt.date(self.year, 1, 1).weekday()
        return ((np.arange(self.n_days) + start) % 7).astype(np.int64)


def build_calendar(year: int) -> OpCalendar:
    """Construct the operating calendar for ``year``."""
    start = dt.date(year, 1, 1)
    n_days = (dt.date(year + 1, 1, 1) - start).days
    dates = [start + dt.timedelta(days=int(i)) for i in range(n_days)]
    dow = np.array([d.weekday() for d in dates], dtype=np.int64)
    weekend = dow >= 5
    us = holidays.country_holidays("US", years=year, observed=True)
    tx = holidays.country_holidays("US", subdiv="TX", years=year, observed=True)
    holiday = np.array([d in us for d in dates], dtype=bool)
    protected = np.zeros(n_days, dtype=bool)
    for i, d in enumerate(dates):
        names: set[str] = set()
        if d in us:
            names.update(n.strip() for n in us.get(d, "").split(";"))
        if d in tx:
            names.update(n.strip() for n in tx.get(d, "").split(";"))
        for n in names:
            base = n.replace(" (observed)", "")
            if base in PEA_RESTRICTED_HOLIDAYS:
                protected[i] = True
    business = ~weekend & ~holiday
    memorial = _holiday_date(us, "Memorial Day")
    labor = _holiday_date(us, "Labor Day")
    summer = np.zeros(n_days, dtype=bool)
    if memorial is not None and labor is not None:
        i0 = (memorial - start).days
        i1 = (labor - start).days
        summer[i0 : i1 + 1] = True
    return OpCalendar(
        year=year,
        n_days=n_days,
        business=business,
        weekend=weekend,
        holiday=holiday,
        summer=summer,
        protected_holiday=protected,
    )


def _holiday_date(cal: holidays.HolidayBase, name: str) -> dt.date | None:
    for d, n in sorted(cal.items()):
        if name in n:
            return d
    return None


def is_night(t: FloatArray) -> BoolArray:
    """Night indicator (22:00-07:00) for times in days."""
    hour = (t - np.floor(t)) * 24.0
    return (hour >= NIGHT_START_H) | (hour < NIGHT_END_H)


def draw_time_of_day(rng: np.random.Generator, night: BoolArray) -> FloatArray:
    """Draw a time of day (hours) uniformly within the requested period."""
    u = rng.random(night.shape[0])
    day_h = NIGHT_END_H + u * DAY_HOURS
    night_h = (NIGHT_START_H + u * NIGHT_HOURS) % 24.0
    return np.where(night, night_h, day_h)


def business_time_to_calendar(cal: OpCalendar, b: FloatArray) -> IntArray:
    """Map positions on the business-day axis to calendar day indices.

    The business-day axis concatenates business days end to end so that
    ``b`` in ``[k, k+1)`` falls on the k-th business day of the year.
    """
    bd = cal.business_days
    k = np.clip(np.floor(b).astype(np.int64), 0, bd.shape[0] - 1)
    return bd[k]


def daylight_hours_mask(cal: OpCalendar, resolution_h: float = 1.0) -> BoolArray:
    """Hourly mask of daylight (07:00-22:00) hours over the year."""
    n = int(round(cal.n_days * 24 / resolution_h))
    hours = (np.arange(n) * resolution_h) % 24.0
    return (hours >= NIGHT_END_H) & (hours < NIGHT_START_H)
