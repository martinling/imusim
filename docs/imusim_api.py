import os.path
from docutils import nodes
import inspect
import imusim
import imusim.all
import functools


def imusim_api_role(role, rawtext, text, lineno, inliner, options={}, content=[]):
    """
    Sphinx extension to link into to IMUSim API.
    """
    uribase = "../api/"

    try:
        # Look for object in imusim.all.
        obj = functools.reduce(lambda x,y: getattr(x,y), [imusim.all] + text.split('.'))

        if inspect.ismodule(obj):
            file = '%s-module.html' % obj.__name__
        elif inspect.isclass(obj):
            file = '%s.%s-class.html' % (obj.__module__, obj.__name__)
        elif inspect.isfunction(obj):
            file = '%s-module.html#%s' % (obj.__module__, obj.__name__)
        elif inspect.ismethod(obj):
            cls = obj.im_class
            file = '%s.%s-class.html#%s' \
                    % (cls.__module__, cls.__name__, obj.__name__)
        elif inspect.isbuiltin(obj):
            if hasattr(obj, '__module__'):
                # Native function
                file = obj.__module__ + '-module.html#' + obj.__name__
            elif hasattr(obj, '__objclass__'):
                # Native method
                cls = obj.__objclass__
                file = '%s.%s-class.html#%s' \
                    % (cls.__module__, cls.__name__, obj.__name__)
            else:
                raise TypeError("Don't know how to document native object " + repr(obj))
        else:
            raise TypeError("Don't know how to document Python object " + repr(obj))
    except AttributeError:
        # Look for object as an imusim submodule.
        __import__("imusim.%s" % text)
        obj = functools.reduce(lambda x,y: getattr(x,y), [imusim] + text.split('.'))
        file = 'imusim.%s-module.html' % text
    except ImportError:
        raise KeyError("Could not find an IMUSim object called '%s'" % text)

    if inspect.ismethod(obj) \
            or (inspect.isbuiltin(obj) and hasattr(obj, '__objclass__')):
        name = obj.__name__
    else:
        name = text

    uri = uribase + file

    node = nodes.reference(rawtext, name, refuri=uri, **options)

    return [node], []


def setup(app):
    app.add_role('api', imusim_api_role)
