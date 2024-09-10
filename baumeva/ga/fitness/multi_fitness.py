from typing import List, Union, Callable, Any
from .base_fitness import BaseFitness
from baumeva.ga import BasePenalty
from baumeva.ga.multi_ga_data import MultiGaData


class MultiFitness(BaseFitness):
    """
    Class for calculating fitness values of one population in case of multiobjective optimization with VEGA.
    Inherits from BaseFitness
    attribute: __idx_opt_value: list of indices of optimization values from object functions.
    """
    __idx_opt_value: List[int] = None

    def __init__(self, obj_function: Union[Callable[[any, List[Union[int, float]]], Union[int, float, tuple]], Callable[
                               [List[Union[int, float]]], Union[int, float, tuple]]],
                 obj_value: Union[int, float, List[Union[int, float]]] = None,
                 input_data: Any = None,
                 penalty: BasePenalty = None,
                 conditions: list = None) -> None:
        """
        Initialize the BaseFitness instance.

        :param obj_function: objective function to be solved, not fitness function! Input argument have to type list.
                             Example:
                                def my_func(x:list):
                                    return x[0]*2 - 1
        :param obj_value: desired objective function value.
        :param input_data: additional information for calculating the value of the objective function.
        :param penalty: subclass of BasePenalty(), initialization before initialization subclass of BaseFitness(),
                               used for conditional optimization. Default: None.
        :param conditions: list of strings (optimizer and conditionals) 3 value can be use: 'optimize', '<=', '!='.
                           Default: None.
                           Example: There is objective function: my_obj_func(x1, x2):
                                                                    return x1**2 + x2**2, 1-x1+x2, x1+x2
                                    my_obj_func returns 3 values, first value to optimize, second value must be <= 0,
                                    third value != 0, so have conditions = ['optimize', '<=', '!=']
                                    dp = DynamicPenalty()
                                    HyperbolaFitness(obj_function=my_func, obj_value=0, penalty=dp,
                                                     conditions=['optimize', '<=', '!='])
        :return: None
        """

        self.num_objectives = conditions.count('optimize')
        super().__init__(obj_function, obj_value, input_data, penalty, conditions)

    def get_fitness_score(self, obj_score: Union[int, float], penalty_value: Union[int, float] = 0) ->\
            Union[int, float]:
        """
        Method for calculating fitness score of individual, will be implemented in child classes.

        :param obj_score: object values of an individual.
        :param penalty_value: fitness values (with respect to optimization values) of an individual.
        :return: fitness scores of the individual.
        """
        pass

    def set_opt_value(self) -> None:
        """
        Method for setting indices of optimization values.
        :return: None.
        """
        self.__idx_opt_value = [i for i in reversed(range(len(self.conditions)))
                                if self.conditions[i] == 'optimize']
        self.conditions = list(filter(lambda a: a != 'optimize', self.conditions))

    def check_task(self) -> None:
        """
        Method for definition type of task: conditional or not.
        :return: None.
        """
        super().check_task()

    def check_input(self, values: List[Union[int, float]]) -> List[Union[int, float]]:
        """
        Checks whether the number of values returned from objective function is correct.
        :param values: values of objective function.
        :return: values.
        """
        return list(values)

    def set_obj_score(self, values: List[Union[int, float]], individ: dict) -> None:
        """
        Gets next objective scores from values of objective function.

        :param values: values of objective function
        :param individ: specimen which gets the objective value.
        :return: None.
        """
        individ['obj_score'] = []
        for idx in self.__idx_opt_value:
            individ['obj_score'].insert(0, values.pop(idx))

    def dominated(self, x: dict, y: dict) -> Union[dict, None]:
        """
        Computes the dominated individual of 2

        :param x: dict, first individual data.
        :param y: dict, second individual data.
        :return: dict, containing the dominated individual of x and y if there is one,
        None if there are no dominated individuals
        """
        x_score = x['obj_score']
        y_score = y['obj_score']

        if x_score == y_score:
            return None

        x_is_dominated = 0

        for i in range(len(y_score)):
            if x_score[i] < y_score[i]:
                if x_is_dominated == 0:
                    x_is_dominated = -1
                elif x_is_dominated == 1:
                    x_is_dominated = 0
                    break
            elif y_score[i] < x_score[i]:
                if x_is_dominated == 0:
                    x_is_dominated = 1
                elif x_is_dominated == -1:
                    x_is_dominated = 0
                    break

        if x_is_dominated == 1:
            return x
        elif x_is_dominated == -1:
            return y
        else:
            return None

    def assign_ranks(self, ga_data: MultiGaData) -> None:
        """
        Assigns a rank to every individual in population equal to 1 + number of individuals dominating it

        :return: None
        """
        for individ in ga_data.population:
            individ['rank'] = 1

        for i in range(len(ga_data.population)):
            for j in range(i+1, len(ga_data.population)):
                inf = self.dominated(ga_data.population[i], ga_data.population[j])
                if inf is not None:
                    inf['rank'] += 1

    def execute(self, ga_data: MultiGaData) -> None:
        """
        Calculate and assign fitness scores to individuals in the population.

        :param ga_data: GaData instance containing population and related data.
        :return: None
        """

        if ga_data.population.is_phenotype:
            ga_data.population.get_phenotype()
            ga_data.population.swap()

        for individ in ga_data.population:
            values = self.calc_obj_func(genotype=individ['genotype'])
            if self._BaseFitness__is_conditional_opt:
                self.set_obj_score(values, individ)
                for i_v, v in enumerate(values):
                    if self.conditions[i_v] == '<=':
                        if v > 0:
                            individ['feasible'] = False
                    else:
                        if v != 0:
                            individ['feasible'] = False

                penalty_value = self.get_penalty_value(values=values, idx_generation=ga_data.idx_generation,
                                                       best_individ=ga_data.best_solution)
                for idx in range(len(individ['obj_score'])):
                    individ['obj_score'][idx] += penalty_value
            else:
                individ['obj_score'] = values

        self.assign_ranks(ga_data)

        for individ in ga_data.population:
            individ['score'] = self.get_fitness_score(individ)

        if ga_data.population.is_phenotype:
            ga_data.population.swap()

        # if ga_data.population.is_phenotype:
        #     ga_data.population.get_phenotype()
        # super().execute(ga_data)