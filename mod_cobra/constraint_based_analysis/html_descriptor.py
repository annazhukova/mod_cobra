from jinja2 import Environment, PackageLoader

__author__ = 'anna'


def describe(template, **kwargs):
    env = Environment(
        loader=PackageLoader('mod_cobra.constraint_based_analysis.cobra_constraint_based_analysis', 'templates'))
    template = env.get_template(template)
    return template.render(**kwargs)
