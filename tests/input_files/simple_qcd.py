################################################################################
#
# Copyright (c) 2009 The MadGraph5_aMC@NLO Development team and Contributors
#
# This file is a part of the MadGraph5_aMC@NLO project, an application which 
# automatically generates Feynman diagrams and matrix elements for arbitrary
# high-energy processes in the Standard Model and beyond.
#
# It is subject to the MadGraph5_aMC@NLO license which should accompany this 
# distribution.
#
# For more information, visit madgraph.phys.ucl.ac.be and amcatnlo.web.cern.ch
#
################################################################################
"""A simple QCD model for unit tests."""

import madgraph.core.base_objects as base_objects
import madgraph.core.color_algebra as color

import copy
# import os

# pjoin = os.path.join

# Particles
#===============================================================================

# A gluon
gluon = base_objects.Particle({
    'name': 'g',
    'antiname': 'g',
    'spin': 3,
    'color': 8,
    'mass': 'zero',
    'width': 'zero',
    'texname': 'g',
    'antitexname': 'g',
    'line': 'curly',
    'charge': 0.,
    'pdg_code': 21,
    'propagating': True,
    'is_part': True,
    'self_antipart': True
})

# Up quark
up = base_objects.Particle({
    'name': 'u',
    'antiname': 'u~',
    'spin': 2,
    'color': 3,
    'mass': 'zero',
    'width': 'zero',
    'texname': 'u',
    'antitexname': '\bar u',
    'line': 'straight',
    'charge': 2. / 3.,
    'pdg_code': 2,
    'propagating': True,
    'is_part': True,
    'self_antipart': False
})
antiup = copy.copy(up)
antiup.set('is_part', False)

# Down quark
down = base_objects.Particle({
    'name': 'd',
    'antiname': 'd~',
    'spin': 2,
    'color': 3,
    'mass': 'zero',
    'width': 'zero',
    'texname': 'u',
    'antitexname': '\bar u',
    'line': 'straight',
    'charge': -1. / 3.,
    'pdg_code': 1,
    'propagating': True,
    'is_part': True,
    'self_antipart': False
})
antidown = copy.copy(down)
antidown.set('is_part', False)

# Photon
photon = base_objects.Particle({
    'name': 'a',
    'antiname': 'a',
    'spin': 3,
    'color': 1,
    'mass': 'zero',
    'width': 'zero',
    'texname': '\gamma',
    'antitexname': '\gamma',
    'line': 'wavy',
    'charge': 0.,
    'pdg_code': 22,
    'propagating': True,
    'is_part': True,
    'self_antipart': True
})

# Higgs
higgs = base_objects.Particle({
    'name': 'h',
    'antiname': 'h',
    'spin': 1,
    'color': 1,
    'mass': 'mh',
    'width': 'wh',
    'texname': 'h',
    'antitexname': 'h',
    'line': 'dashed',
    'charge': 0.,
    'pdg_code': 25,
    'propagating': True,
    'is_part': True,
    'self_antipart': True
})

particles = base_objects.ParticleList([gluon, up, down, photon, higgs])

# Interactions
#===============================================================================

interactions = base_objects.InteractionList()

# g g g
interactions.append(base_objects.Interaction({
    'id': 1,
    'particles': base_objects.ParticleList([gluon, ] * 3),
    'color': [color.ColorString([color.f(0, 1, 2)])],
    'lorentz': ['L1'],
    'couplings': {(0, 0): 'G'},
    'orders': {'QCD': 1}
}))

# g g g g
interactions.append(base_objects.Interaction({
    'id': 2,
    'particles': base_objects.ParticleList([gluon, ] * 4),
    'color': [
        color.ColorString([color.f(-1, 0, 2), color.f(-1, 1, 3)]),
        color.ColorString([color.f(-1, 0, 3), color.f(-1, 1, 2)]),
        color.ColorString([color.f(-1, 0, 1), color.f(-1, 2, 3)])
    ],
    'lorentz': ['L(p1,p2,p3)', 'L(p2,p3,p1)', 'L3'],
    'couplings': {(0, 0): 'G^2', (1, 1): 'G^2', (2, 2): 'G^2'},
    'orders': {'QCD': 2}
}))

# u u~ g
interactions.append(base_objects.Interaction({
    'id': 3,
    'particles': base_objects.ParticleList([
        up, antiup, gluon
    ]),
    'color': [color.ColorString([color.T(2, 0, 1)])],
    'lorentz': ['L1'],
    'couplings': {(0, 0): 'GQQ'},
    'orders': {'QCD': 1}
}))

# d d~ g
interactions.append(base_objects.Interaction({
    'id': 4,
    'particles': base_objects.ParticleList([
        down, antidown, gluon
    ]),
    'color': [color.ColorString([color.T(2, 0, 1)])],
    'lorentz': ['L1'],
    'couplings': {(0, 0): 'GQQ'},
    'orders': {'QCD': 1}
}))

# u u~ y
interactions.append(base_objects.Interaction({
    'id': 5,
    'particles': base_objects.ParticleList([
        up, antiup, photon
    ]),
    'color': [color.ColorString([color.T(0, 1)])],
    'lorentz': ['L1'],
    'couplings': {(0, 0): 'GQED'},
    'orders': {'QED': 1}
}))

# Model
#===============================================================================

model = base_objects.Model({
    'particles': particles,
    'name': "simple_QCD",
    'interactions': interactions
})
