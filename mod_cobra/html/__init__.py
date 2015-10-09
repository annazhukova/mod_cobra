from jinja2 import Environment, PackageLoader

__author__ = 'anna'


def describe(template, **kwargs):
    env = Environment(
        loader=PackageLoader('mod_cobra.html', 'templates'))
    template = env.get_template(template)
    kwargs = {key: (value.decode('utf-8') if isinstance(value, str) else value) for (key, value) in kwargs.iteritems()}
    return template.render(**kwargs)
