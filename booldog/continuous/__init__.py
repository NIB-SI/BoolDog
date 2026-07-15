'''Continuous / semi-quantitative relaxations of Boolean networks.

Exposes :py:class:`~booldog.continuous.semi_quantitative.ContinuousMixin`,
the mixin that gives :py:class:`booldog.BoolDogModel` its
`transform_bool_to_continuous`/`continuous_simulation` methods. The
actual ODE systems
(:py:class:`~booldog.continuous.ode_factory.BooleCubeODE`,
:py:class:`~booldog.continuous.ode_factory.SquadODE`) are built by
:py:mod:`booldog.continuous.ode_factory`.
'''

from .semi_quantitative import ContinuousMixin
