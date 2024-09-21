import datetime
import math
from collections.abc import Sequence
from pathlib import Path
from typing import Literal, Self

from ortools.math_opt.python import mathopt
from pydantic import BaseModel, ConfigDict
from pydanticio import read_records_from_file


class CoordinateRecord(BaseModel):
    model_config = ConfigDict(frozen=True)

    x: float
    y: float


type DistanceType = Literal["euclidian", "manhattan"]


def euclidian_distance(ci: CoordinateRecord, cj: CoordinateRecord) -> float:
    return math.sqrt(pow(ci.x - cj.x, 2) + pow(ci.y - cj.y, 2))


def manhattan_distance(ci: CoordinateRecord, cj: CoordinateRecord) -> float:
    return abs(ci.x - cj.x) + abs(ci.y - cj.y)


class DistanceCalculator:
    def __init__(self, mat: list[list[float]]) -> None:
        self.mat = mat

    def __call__(self, i: int, j: int) -> float:
        return self.mat[i][j]

    def size(self) -> int:
        return len(self.mat)

    @classmethod
    def from_coordinates(
        cls, coords: Sequence[CoordinateRecord], distance_type: DistanceType
    ) -> Self:
        match distance_type:
            case "euclidian":
                dist_func = euclidian_distance
            case "manhattan":
                dist_func = manhattan_distance
            case _:
                raise ValueError(f"Unknown distance type: {distance_type}")
        mat = [[0.0] * len(coords) for _ in range(len(coords))]
        for i, ci in enumerate(coords):
            for j, cj in enumerate(coords):
                mat[i][j] = dist_func(ci, cj)
        return cls(mat)

    @classmethod
    def try_from_coordinate_file(
        cls, path: str | Path, distance_type: DistanceType
    ) -> Self:
        coordinates = read_records_from_file(path, CoordinateRecord)
        return cls.from_coordinates(coordinates, distance_type)


type SolverType = Literal["scip", "highs", "cpsat"]


def solve_tsp(
    dc: DistanceCalculator,
    solver_name: SolverType,
    time_limit: int,
    verbose: bool,
):
    n = dc.size()
    model = mathopt.Model()
    x = [[model.add_binary_variable() for _ in range(n)] for _ in range(n)]
    u = [model.add_integer_variable(lb=0, ub=n - 1) for _ in range(n)]

    obj = sum(dc(i, j) * x[i][j] for i in range(n) for j in range(n))
    model.minimize_linear_objective(obj)

    constraints = []

    # edge constraints
    for i in range(n):
        constraints.append(x[i][i] == 0)
        constraints.append(sum(x[i][j] for j in range(n)) == 1)
        constraints.append(sum(x[j][i] for j in range(n)) == 1)

    # start and end at 0
    constraints.append(u[0] == 0)
    for i in range(1, n):
        constraints.append(u[i] >= 1)

    # subtour elimination by MTZ formulation
    for i in range(1, n):
        for j in range(1, n):
            if i != j:
                constraints.append(u[i] - u[j] <= n * (1 - x[i][j]) - 1)

    # add all constraints to the model
    for c in constraints:
        model.add_linear_constraint(c)

    params = mathopt.SolveParameters(
        time_limit=datetime.timedelta(seconds=time_limit), enable_output=verbose
    )

    match solver_name:
        case "highs":
            solver = mathopt.SolverType.HIGHS
        case "scip":
            solver = mathopt.SolverType.GSCIP
        case "cpsat":
            solver = mathopt.SolverType.CP_SAT
        case _:
            raise ValueError(f"Unknown solver: {solver_name}")

    result = mathopt.solve(model, solver, params=params)

    assert result.termination.reason in (
        mathopt.TerminationReason.OPTIMAL,
        mathopt.TerminationReason.FEASIBLE,
    )

    x_vals = [
        [round(result.variable_values(x[i][j])) for j in range(n)] for i in range(n)
    ]

    for i in range(n):
        assert sum(x_vals[i][j] for j in range(n)) == 1
        assert sum(x_vals[j][i] for j in range(n)) == 1

    tour = [0]
    current = 0
    while True:
        j = x_vals[current].index(1)
        if j == 0:
            break
        tour.append(j)
        current = j

    total_distance = int(result.objective_value())

    return total_distance, tour


def main(
    coordinate_file: Path,
    /,
    distance_type: DistanceType = "euclidian",
    solver: SolverType = "cpsat",
    time_limit: int = 30,
    verbose: bool = False,
):
    dc = DistanceCalculator.try_from_coordinate_file(coordinate_file, distance_type)
    total_distance, tour = solve_tsp(dc, solver, time_limit, verbose)
    print(f"Total distance: {total_distance}")
    print(f"Tour: {tour}")


if __name__ == "__main__":
    import tyro

    tyro.cli(main)
