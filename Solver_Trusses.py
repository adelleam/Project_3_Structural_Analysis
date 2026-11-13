#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 14 14:34:19 2021

@author: kendrick
"""

import numpy as np

# compute unknown displacements 
def ComputeDisplacements(K, F, n_unknowns):
    # extract submatrix of unknowns
    K11 = K[0:n_unknowns,0:n_unknowns]
    F1 = F[0:n_unknowns]
    
    d = np.linalg.solve(K11,F1)
    
    return d

# postprocess the forces at known displacement nodes
def PostprocessReactions(K, d, F, n_unknowns, nodes):
    # These are computed net forces and do not
    # take into account external loads applied
    # at these nodes
    F = np.matmul(K[n_unknowns:,0:n_unknowns], d)
    
    # Postprocess the reactions
    for node in nodes:
        if node.xidx >= n_unknowns:
            node.AddReactionXForce(F[node.xidx-n_unknowns][0] - node.xforce_external)
        if node.yidx >= n_unknowns:
            node.AddReactionYForce(F[node.yidx-n_unknowns][0] - node.yforce_external)
        
    return F

# determine internal member loads
def ComputeMemberForces(bars):
    # COMPLETE THIS FUNCTION
    # Compute member forces for all bars using equation 14-23 
    for bar in bars:
        #bar properties
        E = bar.E
        A = bar.A
        L = bar.Length()
        
        #direction of bar
        lambdax, lambday = bar.LambdaTerms()
        
        #horizontal/vertical movement of starting node
        u1x = bar.init_node.xdisp
        u1y = bar.init_node.ydisp
        
        #get movements of ending node
        u2x = bar.end_node.xdisp
        u2y = bar.end_node.ydisp
        
        #calculate internal force
        axial_force = (E*A/L) * (-lambdax*u1x - lambday*u1y + lambdax*u2x + lambday*u2y)
        
        bar.axial_load = axial_force
    
# compute the normal stresses
def ComputeNormalStresses(bars):
    # COMPLETE THIS FUNCTION
    # Compute normal stress for all bars
    
    for bar in bars:
        
        #axial force in bar
        axial_force = bar.axial_load
        
        #cross-sectional area of bar
        A = bar.A
        
        #calculate normal stress (F/A)
        normal_stress = axial_force / A
        
        bar.normal_stress = normal_stress


# compute the critical buckling load of a member
def ComputeBucklingLoad(bars):
    # COMPLETE THIS FUNCTION
    # Compute critical buckling load for all bars
    
    K = 1
    
    for bar in bars:
        
        #material stiffness and length
        E = bar.E
        L = bar.Length()
        
        #strong and weak axis moments of inertia
        Istrong = bar.GetStrongSecondMoment()
        Iweak = bar.GetWeakSecondMoment()
        
        #calculate buckling with strong and weak
        Pstrong = (np.pi**2 * E * Istrong) / (K * L)**2
        Pweak = (np.pi**2 * E * Iweak) / (K*L)**2
        
        #Choose smaller load
        bar.buckling_load = min(Pstrong, Pweak)
