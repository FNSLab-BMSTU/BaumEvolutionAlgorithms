from .multi_ga_data import MultiGaData
from .new_generation import NewGeneration


class MultiNewGeneration(NewGeneration):
    """
    Class for creating a new generation of individuals in a genetic algorithm in case of multiobjective optimization.
    """
    def __init__(self, transfer_parents: str = 'best', num_elites: int = 1) -> None:
        """
        Initialize the MultiNewGeneration instance.

        :param transfer_parents: strategy for transferring parents to the next generation.
        :param num_elites: number of elites in case of the 'best' strategy.
        :return: None
        """
        self.num_elites = num_elites
        super().__init__(transfer_parents)

    def add_best(self, ga_data: MultiGaData, num_elites):
        """
        Add best parent individuals to the offspring (elitism strategy).

        :param ga_data: MultiGaData instance containing population and related data.
        :param num_elites: number of parent individuals to add.
        :return: list of elites to add to the population.
        """
        elites = ga_data.population[:num_elites]

        return elites

    def execute(self, ga_data: MultiGaData) -> None:
        """
        Execute the new generation creation process.

        :param ga_data: MultiGaData instance containing population and related data.
        :return: None
        """
        if ga_data.population.is_phenotype:
            ga_data.population.get_phenotype()

        super().execute(ga_data)
