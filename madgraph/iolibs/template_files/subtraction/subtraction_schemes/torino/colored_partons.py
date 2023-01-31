"""This function returns the colored partons of a reduced_process as legs"""
import madgraph.core.subtraction as sub

def get_colored_legs(reduced_process):
    """ returns a dictionary {leg_number : leg}
        with only colored legs in it"""

    c_legs_numbers = get_colored_leg_numbers(reduced_process).keys()
    c_legs = {}
    for l in reduced_process.get('legs'):
        if l['number'] in c_legs_numbers:
            c_legs[l['number']] = sub.SubtractionLeg(l)

    return c_legs

def get_colored_leg_numbers(reduced_process):
    """ Returns a dictionary of the form:
        { leg_number : ( SU3_representation, leg_state)  }
        where leg_state is either SubtractionLeg.INITIAL or SubtractionLeg.FINAL and SU3_representation
        is 3 for a triplet, -3 for an antitriplet, 8 for an octet... etc. (legs not present in this dictionary are
        colorless by construction).
    """

# Build a dictionary of the form {leg_number : (color_representation, state [in/outgoing])}
    colored_parton_numbers = {}
    model = reduced_process.get('model')
    for leg in reduced_process.get('legs'):
        leg_color_quantum_number = model.get_particle(leg.get('id')).get('color')
        if leg_color_quantum_number==1:
            continue
        colored_parton_numbers[leg.get('number')] = (leg_color_quantum_number, leg.get('state'))

    return colored_parton_numbers