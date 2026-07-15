'''Result objects returned by simulations (not by ``BoolDogModel`` itself):
:class:`~booldog.simulation_result.continuous_result.ContinuousSimulationResult`
from :py:meth:`~booldog.continuous.semi_quantitative.ContinuousMixin.continuous_simulation`,
and :class:`~booldog.simulation_result.boolean_result.BooleanSimulationResult` /
:class:`~booldog.simulation_result.boolean_result.BooleanStateSpace` from
:py:meth:`~booldog.boolean.boolean.BooleanNetworkMixin.boolean_simulation`.
'''

from .continuous_result import ContinuousSimulationResult
from .boolean_result import BooleanSimulationResult, BooleanStateSpace
