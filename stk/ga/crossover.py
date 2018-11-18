"""
Defines crossover operations via the :class:`Crossover` class.

.. _`adding crossover functions`:

Extending stk: Adding crossover functions.
------------------------------------------

If a new crossover operation is to be added to ``stk`` it should be
added as a method in the :class:`Crossover` class defined in this
module. The only requirement is that the function returns the offpsring
in a :class:`.GAPopulation` instance. The function will be called with
each parent selected by the crossover selection function in a
separate argument. For example if you wish to define a function
which carries out crossover across two parents

.. code-block:: python

    def crossover_fn_example(self, macro_mol1, macro_mol2, param1):
        ...

or with a variable number of parents

.. code-block:: python

    def crossover_fn_example2(self, *macro_mol, param2):
        ...

and of course any where in between is valid too. Just make sure that
the crossover selection function selects the number of parents that the
crossover function expects to receive when running the GA.

If the crossover function does not fit neatly into a single function
make sure that any helper functions are private, i.e. that their names
start with a leading underscore.

"""

import logging
from collections import Counter, defaultdict
import numpy as np
from itertools import islice, product

from .ga_population import GAPopulation
from .plotting import plot_counter
from ..utilities import dedupe


logger = logging.getLogger(__name__)


class Crossover:
    """
    Carries out crossover operations on the population.

    Instances of :class:`.GAPopulation` delegate crossover operations
    to instances of this class. They do this by calling

    >>> offspring_pop = pop.gen_offspring()

    where ``offspring_pop`` is a new :class:`.GAPopulation`, holding
    molecules generated by performing crossover operations on members
    of ``pop``. This class uses the :class:`.Selection` instance in
    ``pop.ga_tools.selection`` to select parents used in crossover.

    Attributes
    ----------
    funcs : :class:`list` of :class:`.FunctionData`
        This lists holds all the crossover functions which are to be
        used. One will be chosen at random when a crossover operation
        is to be performed. The likelihood that each is selected is
        given by :attr:`weights`.

    num_crossovers : :class:`int`
        The number of crossover operations performed each time
        :meth:`.GAPopulation.gen_offspring` is called.

    weights : :class:`list` of :class:`float`
        Each float corresponds to the probability of selecting the
        crossover function in :attr:`funcs` at the corresponding index.
        For example,

        .. code-block:: python

            selection = Selection(funcs=[FunctionData('one'),
                                         FunctionData('two')],
                                  num_crossovers=3,
                                  weights=[0.3, 0.7])

        means that the crossover function called "one" has a
        probability of ``0.3`` of being used, while the crossover
        function called "two" has a probability of ``0.7`` of being
        used.

        This means entries in this list must sum to 1 and the number of
        entries must be the same as in :attr:`funcs`. Defaults to
        ``None``, which means all crossover functions have an equal
        probability of selection.

    """

    def __init__(self, funcs, num_crossovers, weights=None):
        """
        Intializes a :class:`Crossover` object.

        Parameters
        ----------
        funcs : :class:`list` of :class:`.FunctionData`
            This lists holds all the crossover functions which are to
            be used. One will be chosen at random when a crossover
            operation is to be performed. The likelihood that each is
            selected is given by :attr:`weights`.

        num_crossovers : :class:`int`
            The number of crossover operations performed each time
            :meth:`.GAPopulation.gen_offspring` is called.

        weights : :class:`list` of :class:`float`, optional
            Each float corresponds to the probability of selecting the
            crossover function in :attr:`funcs` at the corresponding
            index. This means entries in this list must sum to 1 and
            the number of entries must be the same as in :attr:`funcs`.
            Defaults to ``None``, which means all crossover functions
            have an equal probability of selection.

        """

        self.funcs = funcs
        self.weights = weights
        self.num_crossovers = num_crossovers

    def __call__(self, population, counter_path=''):
        """
        Carries out crossover operations on `population`.

        This function selects members of `population` and crosses
        them until either all possible parents have been crossed or the
        required number of successful crossover operations has been
        performed.

        The offspring generated are returned together in a
        :class:`.GAPopulation` instance. Any molecules that are created
        via crossover and match a molecule present in the original
        population are removed.

        Parameters
        ----------
        population : :class:`.GAPopulation`
            The population instance who's members are to be crossed.

        counter_path : :class:`str`, optional
            The name of the ``.png`` file showing which members were
            selected for crossover. If ``''``, then no file is made.

        Returns
        -------
        :class:`.GAPopulation`
            A population with all the offspring generated held in its
            :attr:`~.Population.members` attribute. This does not
            include offspring which correspond to molecules already
            present in `population`.

        """

        offspring_pop = GAPopulation(ga_tools=population.ga_tools)
        counter = Counter()

        parent_pool = islice(population.select('crossover'),
                             self.num_crossovers)
        for i, parents in enumerate(parent_pool, 1):
            logger.info('Crossover number {}. Finish when {}.'.format(
                                           i, self.num_crossovers))
            counter.update(parents)
            # Get the crossover function.
            func_data = np.random.choice(self.funcs, p=self.weights)
            func = getattr(self, func_data.name)
            logger.info(f'Using {func.__name__}.')

            try:
                # Apply the crossover function and supply any
                # additional arguments to it.
                offspring = func(*parents, **func_data.params)

                # Print the names of offspring which have been returned
                # from the cache.
                for o in offspring:
                    if o.name:
                        logger.debug(('Offspring "{}" retrieved '
                                      'from cache.').format(o.name))

                # Add the new offspring to the offspring population.
                offspring_pop.add_members(offspring)

            except Exception:
                errormsg = ('Crossover function "{}()" failed on '
                            'molecules PARENTS.').format(
                            func_data.name)

                pnames = ' and '.join('"{}"'.format(p.name) for
                                      p in parents)
                errormsg = errormsg.replace('PARENTS', pnames)
                logger.error(errormsg, exc_info=True)

        # Make sure that only original molecules are left in the
        # offspring population.
        offspring_pop -= population

        if counter_path:
            # Update counter with unselected members and plot counter.
            for member in population:
                if member not in counter.keys():
                    counter.update({member: 0})
            plot_counter(counter, counter_path)

        return offspring_pop

    def genetic_recombination(self, *macro_mols, key, n_offspring=1):
        """
        Recombine building blocks using biological systems as a model.

        Overall, this function mimics how animals and plants inherit
        DNA from their parents, except generalized to work with any
        number of parents. First it is worth discussing some
        terminology. A gene is a the smallest packet of genetic
        information. In animals, each gene can have multiple alleles.
        For example, there is a gene for hair color, and individual
        alleles for black, red, brown, etc. hair. This means that every
        person has a gene for hair color, but a person with black hair
        will have the black hair allele and a person with red hair will
        have the red hair allele. When two parents produce an
        offspring, the offspring will have a hair color gene and will
        inherit the allele of one of the parents at random. Therefore,
        if you have two parents, one with black hair and one with red
        hair, the offspring will either have black or red hair,
        depending on which allele they inherit.

        In ``stk`` molecules, each building block represents an allele.
        The question is, which gene is each building block an allele
        of? To answer that, let's first construct a couple of
        building block molecules

        .. code-block:: python

            bb1 = StructUnit2('filename1.mol', 'amine')
            bb2 = StructUnit3('filename2.mol', 'aldehyde')
            bb3 = StructUnit2('filename3.mol', 'aldehyde')
            bb4 = StructUnit3('filename4.mol', 'amine')

        We can define a function which analyzes a building block
        molecule and returns the gene it belongs to, for example

        .. code-block:: python

            def determine_gene(building_block):
                return building_block.func_grp.name

        Here, we can see that the gene to which each building block
        molecule belongs is given by the functional group name.
        Therefore there is an ``'amine'`` gene which has two alleles
        ``bb1`` and ``bb4`` and there is an ``'aldehyde'`` gene which
        has two alleles ``bb2`` and ``bb3``.

        Alternatively, we could have defined a function such as

        .. code-block:: python

            def determine_gene(building_block):
                return building_block.__class__.__name__

        Now we can see that we end up with the gene called
        ``'StructUnit2'`` which has two alleles ``bb1`` and ``bb3``
        and a second gene called ``'StructUnit3'`` which has the
        alleles ``bb2`` and ``bb4``.

        To produce offspring molecules, this function categorizes
        each building block of the parent molecules into genes using
        the `key` parameter. Then, to generate a single offspring, it
        picks a random building block for every gene. The picked
        building blocks are used to construct the offspring. The
        topoogy of the offspring one of the parent's topologies, also
        selected at random. For obvious reasons, this approach works
        with any number of parents.

        Parameters
        ----------
        macro_mols : :class:`.MacroMolecule`
            The parent molecules.

        key : :class:`function`
            A function, which takes a :class:`.StructUnit` object
            and returns its gene or category. To produce an offspring,
            one of the building blocks from each category is picked
            at random.

        n_offspring : :class:`int`
            The maximum number of offspring to create.

        Returns
        -------
        :class:`.GAPopulation`
            A population holding all the generated offspring.

        """

        cls = macro_mols[0].__class__
        topologies = [macro_mol.topology for macro_mol in macro_mols]

        genes = defaultdict(set)
        for macro_mol in macro_mols:
            for allele in macro_mol.building_blocks:
                genes[key(allele)].add(allele)

        genes = {gene: np.random.permutation(list(alleles))
                 for gene, alleles in genes.items()}

        offspring_pop = GAPopulation()
        for i, building_blocks in enumerate(product(*genes.values())):
            topology = np.random.choice(topologies)
            offspring = cls(building_blocks, topology)
            offspring_pop.members.append(offspring)

            if i == n_offspring:
                return offspring_pop

        return offspring_pop

    def bb_lk_exchange(self, macro_mol1, macro_mol2):
        """
        Exchanges the building blocks and linkers of cages.

        This operation is basically::

            bb1-lk1 + bb2-lk2 --> bb1-lk2 + bb2-lk1,

        where bb-lk represents a building block - linker combination
        of a cage.

        If the parent cages do not have the same topology, then pairs
        of offspring are created for each topology. This means that
        there may be up to ``4`` offspring.

        Parameters
        ----------
        macro_mol1 : :class:`.MacroMolecule`
            The first parent cage. Its building-block* and linker are
            combined with those of `cage2` to form new cages.

        macro_mol2 : :class:`.MacroMolecule`
            The second parent cage. Its building-block* and linker are
            combined with those of `cage1` to form new cages.

        Returns
        -------
        :class:`.GAPopulation`
            A population of all the offspring generated by crossover of
            `macro_mol1` with `macro_mol2`.

        """

        Cage = macro_mol1.__class__

        # Make a variable for each building-block* and linker of each
        # each cage. Make a set consisting of topologies of the cages
        # provided as arguments - this automatically removes copies.
        # For each topology create two offspring cages by combining the
        # building-block* of one cage with the linker of the other.
        # Place each new cage into a ``GAPopulation`` instance and return
        # that.

        _, c1_lk = max(zip(macro_mol1.bb_counter.values(),
                           macro_mol1.bb_counter.keys()))
        _, c1_bb = min(zip(macro_mol1.bb_counter.values(),
                           macro_mol1.bb_counter.keys()))

        _, c2_lk = max(zip(macro_mol2.bb_counter.values(),
                           macro_mol2.bb_counter.keys()))
        _, c2_bb = min(zip(macro_mol2.bb_counter.values(),
                           macro_mol2.bb_counter.keys()))

        offspring_pop = GAPopulation()
        # For each topology create a new pair of offspring using the
        # building block pairings determined earlier.
        topologies = (x.topology for x in (macro_mol1, macro_mol2))
        for topology in topologies:
            offspring1 = Cage((c1_lk, c2_bb), topology)
            offspring2 = Cage((c2_lk, c1_bb), topology)
            offspring_pop.add_members((offspring1, offspring2))

        return offspring_pop

    def jumble(self,
               macro_mol1,
               macro_mol2,
               n_offspring_building_blocks,
               n_offspring,
               allow_duplicate_building_blocks=False):
        """
        Randomly distributes all building blocks among offspring.

        Puts all the building blocks from each parent into one big pot
        and samples ``N`` building block sets, where each set has a
        size of ``M``. Each sampled set is used to create one new
        offspring. ``N`` is given by `n_offspring` and ``M`` is given
        by `n_offspring_building_blocks`.

        The offspring inherit the topology of one of the parents,
        chosen at random.

        Parameters
        ----------
        macro_mol1 : :class:`.MacroMolecule`
            The first parent.

        macro_mol2 : :class:`.MacroMolecule`
            The second parent.

        n_offspring_building_blocks : :class:`int`
            The number of building blocks each offspring is made from.

        n_offspring : :class:`int`
            The number of offspring to produce.

        allow_duplicate_building_blocks : :class:`bool`, optional
            Indicates whether the building blocks used to construct the
            offspring must all be unique.

        Returns
        -------
        :class:`.GAPopulation`
            A population of all the offspring generated by crossover of
            `macro_mol1` with `macro_mol2`.

        """

        building_blocks = list(dedupe(macro_mol1.building_blocks +
                                      macro_mol2.building_blocks))
        parent_topologies = [macro_mol1.topology, macro_mol2.topology]

        pop = GAPopulation()
        for i in range(n_offspring):
            offspring_topology = np.random.choice(parent_topologies)
            bbs = np.random.choice(
                        building_blocks,
                        size=n_offspring_building_blocks,
                        replace=allow_duplicate_building_blocks)
            offspring = macro_mol1.__class__(bbs.tolist(),
                                             offspring_topology)
            pop.members.append(offspring)
        return pop
