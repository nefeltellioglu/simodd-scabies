from random import Random
from unittest import TestCase

from population.individual import Individual
from population.pop_base import Population


__author__ = 'ngeard'


class TestPopulationIndividual(TestCase):
    def setUp(self):
        self.P = Population(Individual)

    def test_add_individual(self):
        # add an individual to an empty population
        new_age = 10
        new_ind = self.P.add_individual(new_age)
        # check that individual is stored correctly and next_id is incremented
        self.assertEqual(len(self.P.I), 1)
        self.assertEqual(new_ind.ID, 0)
        self.assertEqual(self.P.next_id['individual'], 1)
        self.assertEqual(self.P.I[0].age, new_age)

    def test_remove_individual(self):
        new_age = 10
        new_ind = self.P.add_individual(new_age)
        self.assertEqual(len(self.P.I), 1)
        self.P.remove_individual(new_ind)
        self.assertEqual(len(self.P.I), 0)
        self.assertEqual(self.P.next_id['individual'], 1)

    def test_age_population(self):
        ages = [5, 7, 12, 17, 25]
        # add individuals of varying ages
        for x in ages:
            self.P.add_individual(x)
        ages = [6, 8, 13, 18, 26]
        self.P.age_population(365)
        new_ages = [x.age for x in list(self.P.I.values())]
        self.assertEqual(ages, new_ages)

    def test_individuals_by_age(self):
        ages = [5, 7, 12, 17, 25]
        # add individuals of varying ages
        for x in ages:
            self.P.add_individual(x)
        # ensure function returns appropriate I
        returned_individuals = self.P.individuals_by_age(1, 100)
        returned_ages = [x.age for x in returned_individuals]
        self.assertEqual(len(returned_individuals), 5)
        self.assertEqual(returned_ages, ages)
        returned_individuals = self.P.individuals_by_age(1, 10)
        returned_ages = [x.age for x in returned_individuals]
        self.assertEqual(len(returned_individuals), 2)
        self.assertEqual(returned_ages, [5, 7])
        self.assertEqual(len(self.P.individuals_by_age(1, 3)), 0)


class TestPopulationGroups(TestCase):
    def setUp(self):
        self.P = Population(Individual)
        ages = [5, 7, 12, 17, 25]
        # add I of varying ages
        for x in ages:
            self.P.add_individual(x)

    def test_init_group_type(self):
        label = 'group'
        self.P.init_group_type(label)
        self.assertEqual(len(self.P.groups), 1)
        self.assertTrue(label in self.P.groups)

    def test_add_group(self):
        label = 'group'
        members = [self.P.I[1], self.P.I[3]]
        self.P.init_group_type(label)
        new_id = self.P.add_group(label, members)
        self.assertTrue(new_id in self.P.groups[label])
        self.assertEqual(len(self.P.groups[label]), 1)
        self.assertEqual(new_id, 0)
        self.assertEqual(len(self.P.groups[label][new_id]), 2)
        self.assertEqual(self.P.I[1].groups[label], new_id)

    def test_remove_group(self):
        label = 'group'
        self.P.init_group_type(label)
        members = [self.P.I[1], self.P.I[3]]
        new_id = self.P.add_group(label, members)
        self.P.remove_group(label, new_id)
        self.assertEqual(len(self.P.groups[label]), 0)
        self.assertEqual(len(self.P.I[1].groups), 0)

    def test_groups_by_size(self):
        label = 'group'
        self.P.init_group_type(label)
        id1 = self.P.add_group(label, [self.P.I[0]])
        id2 = self.P.add_group(label, [self.P.I[1], self.P.I[2]])
        id3 = self.P.add_group(label, [self.P.I[3]])
        self.assertEqual(self.P.groups_by_size(label, 1), [id1, id3])
        self.assertEqual(self.P.groups_by_size(label, 2), [id2])
        self.assertEqual(self.P.groups_by_size(label, 3), [])

    def test_individuals_by_group_size(self):
        label = 'group'
        self.P.init_group_type(label)
        self.P.add_group(label, [self.P.I[0]])
        self.P.add_group(label, [self.P.I[1], self.P.I[2]])
        self.P.add_group(label, [self.P.I[3]])
        self.assertEqual([x.ID for x in self.P.individuals_by_group_size(label, 1)], [0, 3])
        self.assertEqual([x.ID for x in self.P.individuals_by_group_size(label, 2)], [1, 2])
        self.assertEqual([x.ID for x in self.P.individuals_by_group_size(label, 3)], [])

    def test_gen_age_structured_pop(self):
        rng = Random(1234)
        age_probs = [[0.1,[1]],[0.3,[2]],[0.4,[3]],[0.2,[4]]]
        # create a population of 10 with given probabilities
        self.P.gen_age_structured_pop(10, age_probs, rng)
        print([x.age for x in self.P.I.values()])
        self.assertEqual([x.age for x in self.P.I.values()], [5, 7, 12, 17, 25, 4, 3, 1, 4, 4, 3, 3, 1, 3, 2])

    def test_allocate_groups_by_age(self):
        rng = Random()
        label = 'group'
        # create five groups of size 1
        self.P.allocate_groups_by_age(label, [[1., [1]]], (1, 100), rng)
        self.assertEqual(len(self.P.groups[label]), 5)
        # remove all groups
        for g in list(self.P.groups[label].keys()):
            self.P.remove_group(label, g)
        # create two groups of size two
        label = 'group2'
        self.P.allocate_groups_by_age(label, [[1., [2]]], (1, 20), rng)
        self.assertEqual(len(self.P.groups[label]), 2)

    def test_add_individuals_to_group(self):
        label = 'group'
        self.P.init_group_type(label)
        group_id = self.P.add_group(label, [self.P.I[0], self.P.I[1], self.P.I[2]])
        self.P.add_individuals_to_group(label, group_id, [self.P.I[3]])
        self.assertEqual(len(self.P.groups[label][group_id]), 4)
        self.assertTrue(label in self.P.I[3].groups)

    def test_remove_individual_from_group(self):
        label = 'group'
        self.P.init_group_type(label)
        self.P.add_group(label, [self.P.I[0], self.P.I[1], self.P.I[2]])
        group_id = self.P.remove_individual_from_group(label, self.P.I[1])
        self.assertEqual(len(self.P.groups[label][group_id]), 2)
        self.assertFalse(label in self.P.I[1].groups)

    def test_remove_individual(self):
        label = 'group'
        self.P.init_group_type(label)
        group_id = self.P.add_group(label, [self.P.I[0], self.P.I[1], self.P.I[2]])
        self.P.remove_individual(self.P.I[0])
        self.assertEqual(len(self.P.I), 4)
        self.assertEqual(len(self.P.groups[label][group_id]), 2)
