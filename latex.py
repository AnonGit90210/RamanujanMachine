"""The following functions are used to generate a latex document and a rendered pdf."""

from pylatex import Document, Section, Alignat

def generate_latex(filename, eqns=[]):
    """Creates a latex document, with a title and the equations.
    eqns - a list of latex math expressions, each of an equation."""
    doc = Document()

    with doc.create(Section('Automatic Conjectures')):
        doc.append('These are the conjectures detected by the algorithm.')

        for eqn in eqns:
            with doc.create(Alignat(numbering=False, escape=False)) as agn:
                agn.append(eqn)

    doc.generate_pdf(filename, clean_tex=False)


def latex_cont_frac(a, b, current_expression=''):
    """Generates a ContFrac latex expression from a, b.
    a - a list of (numerical) values of a for iterative depths.
    b - the same.
    current_expression - for inner use."""
    # The expression is built from bottom (the most inner fraction) up (to the outer fraction).
    # current_expression saves the so far built expression, and is returned at the end
    if current_expression == '':
        current_expression = str(a[-1]) + ' + \dots'

    if len(a) > 1:
        new_iteration = r'{0} + \frac{{ {1} }} {{ {2} }}'.format(a[-2], b[-1], current_expression)
        return latex_cont_frac(a[:-1], b[:-1], new_iteration)
    else:
        return current_expression
    

# if __name__ == '__main__':
#     print(latex_cont_frac([2, 4, 1, 5], [3, 2, 7, 6]))
