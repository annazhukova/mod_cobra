from jinja2 import Environment, PackageLoader

__author__ = 'anna'


def describe(template, **kwargs):
    env = Environment(
        loader=PackageLoader('mod_cobra.html', 'templates'))
    template = env.get_template(template)
    return template.render(**kwargs)
