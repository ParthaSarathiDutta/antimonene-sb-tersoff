"""
Representative evaluator for the antimonene Tersoff fitting workflow.

Users should modify this file according to the properties and objective
functions required for their own fitting problem. 
"""

from dataclasses import dataclass


@dataclass(frozen=True)
class Stage:
    name: str
    initial_tolerance_percent: float
    note: str = ""


# Hierarchical order used in the manuscript.
# No weights are assigned between these five target-property stages.
STAGES = [
    Stage("lattice", 3.0, "Both antimonene polymorphs"),
    Stage("cohesive_energy", 3.0, "Both antimonene polymorphs"),
    Stage("equation_of_state", 30.0, "Both antimonene polymorphs"),
    Stage("elastic_constants", 30.0, "Both antimonene polymorphs"),
    Stage(
        "phonon_dispersion",
        30.0,
        "Beta-antimonene only; acoustic branches were given highest importance",
    ),
]


def evaluate_candidate(params):
    """
    Representative hierarchical evaluator.

    Replace the placeholder property calculations below with the user's own
    simulation/DFT/MD workflow.

    Important details from this study:
      1. A candidate progresses only after satisfying the current-stage
         checkpoint.
      2. The listed tolerances are initial search criteria, not immutable
         thresholds; they may be adjusted between searches.
      3. There are no weights between the five property stages.
      4. Within one property, individual features may contribute differently
         to that property's objective score. For phonons, the acoustic
         branches were prioritized.
    """

    for stage in STAGES:
        predicted = calculate_property(params, stage.name)
        target = load_target(stage.name)

        score = property_objective(stage.name, predicted, target)

        if not passes_checkpoint(
            stage.name,
            predicted,
            target,
            stage.initial_tolerance_percent,
        ):
            return {
                "status": "rejected",
                "failed_stage": stage.name,
                "score": score,
            }

    return {
        "status": "completed",
        "score": final_objective(params),
    }


# -------------------------------------------------------------------------
# Replace these placeholders with problem-specific implementations.
# -------------------------------------------------------------------------

def calculate_property(params, property_name):
    raise NotImplementedError(
        "Implement the property calculation required for your fitting problem."
    )


def load_target(property_name):
    raise NotImplementedError(
        "Load the corresponding reference/training data."
    )


def property_objective(property_name, predicted, target):
    raise NotImplementedError(
        "Define the objective metric for each property. "
        "Sub-features may be weighted within a property if needed."
    )


def passes_checkpoint(property_name, predicted, target, tolerance_percent):
    raise NotImplementedError(
        "Implement the stage-specific checkpoint/error metric."
    )


def final_objective(params):
    raise NotImplementedError(
        "Return the final objective score used to rank completed candidates."
    )
