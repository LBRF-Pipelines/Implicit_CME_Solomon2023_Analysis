__author__ = 'jono'
from klibs.KLIndependentVariable import IndependentVariableSet, IndependentVariable

TraceLab_ind_vars = IndependentVariableSet()

TraceLab_ind_vars.add_variable("animate_time", int)
TraceLab_ind_vars.add_variable("figure_name", str)

TraceLab_ind_vars['animate_time'].add_values(500, 1500, 2500)
#TraceLab_ind_vars['animate_time'].add_values(500, 1500, 2500)
#TraceLab_ind_vars['animate_time'].add_values(10000, 20000)

# 500, 1500, 2500

# this is the default figure set - if you don't name your figure set it'll do this at 50:50
TraceLab_ind_vars['figure_name'].add_values("random", "template_1477090164.31")