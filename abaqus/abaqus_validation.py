from abaqus import *
from abaqusConstants import *
from collections import OrderedDict
backwardCompatibility.setValues(includeDeprecated=True,
                                reportDeprecated=False)

import sketch
import part
import numpy as np
import regionToolset

if __name__ == "__main__":
    modelString = 'bend_ac{0:02d}_at{1:02d}_E{2:04d}_n{3:02d}'
    inpString = 'bend_ac{0:.1f}_at{1:.1f}_L{2:05.2f}_W{3:05.2f}_abq.inp'
    enableFullyPlastic = True
    t = 1.0
    Sys = 1.0
    nu = 0.3
    displacementDict = { 'bend_ac02_at02_E0100_n03': 0.97, # node 37
                         'bend_ac02_at02_E0100_n20': 0.84,
                         'bend_ac02_at02_E0200_n03': 0.97, # node 37
                         'bend_ac02_at02_E0200_n20': 0.84,
                         'bend_ac02_at02_E0500_n03': 0.18,
                         'bend_ac02_at02_E0500_n20': 0.16,
                         'bend_ac02_at02_E1000_n03': 0.09,
                         'bend_ac02_at02_E1000_n20': 0.08,
                         'bend_ac02_at08_E0100_n03': 11.0, # node 18903
                         'bend_ac02_at08_E0100_n20': 9.8,
                         'bend_ac02_at08_E0200_n03': 11.0, # node 18903
                         'bend_ac02_at08_E0200_n20': 9.8,
                         'bend_ac02_at08_E0500_n03': 3.0,
                         'bend_ac02_at08_E0500_n20': 2.6,
                         'bend_ac02_at08_E1000_n03': 1.5,
                         'bend_ac02_at08_E1000_n20': 1.3,
                         'bend_ac10_at02_E0100_n03': 1.0, # node 37
                         'bend_ac10_at02_E0100_n20': 0.95,
                         'bend_ac10_at02_E0200_n03': 1.0, # node 37
                         'bend_ac10_at02_E0200_n20': 0.95,
                         'bend_ac10_at02_E0500_n03': 0.18,
                         'bend_ac10_at02_E0500_n20': 0.16,
                         'bend_ac10_at02_E1000_n03': 0.09,
                         'bend_ac10_at02_E1000_n20': 0.08,
                         'bend_ac10_at08_E0100_n03': 1.3, # node 37
                         'bend_ac10_at08_E0100_n20': 1.1,
                         'bend_ac10_at08_E0200_n03': 1.3, # node 37
                         'bend_ac10_at08_E0200_n20': 1.1,
                         'bend_ac10_at08_E0500_n03': 0.24,
                         'bend_ac10_at08_E0500_n20': 0.24,
                         'bend_ac10_at08_E1000_n03': 0.12,
                         'bend_ac10_at08_E1000_n20': 0.12,
                         }
    for ac in [0.2, 1.0, ]:
        for at in [0.2, 0.8, ]:
            a = at*t
            c = a/ac
            if (5*c > 5*t):
                W = 5*c
            else:
                W = 5*t
            L = 1.1*(2*W)
            S_inner = W
            S_outer = 2*W

            for E in [100, 200, 500, 1000, ]:
                for n in [3, 20]:
                    modelName = modelString.format(int(10*ac), int(10*at),
                                                   E, n)
                    inpName = inpString.format(ac, at, L, W)
                    myModel = mdb.ModelFromInputFile(name=modelName,
                                                     inputFileName=inpName)
                    M = mdb.models[modelName]
                    # Fix material properties
                    if enableFullyPlastic:
                        del M.materials['MAT1'].elastic
                        del M.materials['MAT1'].plastic
                        M.materials['MAT1'].DeformationPlasticity(
                            table=((E, 0.3, 1.0, n, 0.5), ))
                    else:
                        # - Elastic
                        elasticTable = [[E, nu],]
                        M.materials['MAT1'].Elastic(table=elasticTable)
                        # - Plastic
                        plastic_strain = 0.001*np.arange(1, 9)
                        plastic_strain = np.append(plastic_strain,
                                                   plastic_strain[-1]+0.005*np.arange(1, 5))
                        plastic_strain = np.append(plastic_strain,
                                                   plastic_strain[-1]+0.01*np.arange(1, 100))
                        strain_ys = Sys/E
                        strain_powerlaw = strain_ys+plastic_strain
                        stress_powerlaw = Sys*np.power(strain_powerlaw/strain_ys, 1.0/n)
                        plasticTable = [[Sys, 0.0],]
                        for i in range(len(plastic_strain)):
                            plasticTable.append([stress_powerlaw[i], plastic_strain[i]])
                        M.materials['MAT1'].Plastic(table=plasticTable)

                    # Create rigid wall to constrain crack face
                    # - Geometry
                    M.ConstrainedSketch(name='__profile__', sheetSize=200.0)
                    M.sketches['__profile__'].rectangle(point1=(0.0, 0.0),
                                                                point2=(1.0, S_outer))
                    M.Part(dimensionality=THREE_D,
                                               name='Plane',
                                               type=ANALYTIC_RIGID_SURFACE)
                    M.parts['Plane'].AnalyticRigidSurfExtrude(depth=W,
                                                              sketch=M.sketches['__profile__'])
                    M.parts['Plane'].ReferencePoint(point=(0.0, 0.0, 0.0))
                    M.rootAssembly.Instance(dependent=ON, 
                                            name='Plane-1',
                                            part=M.parts['Plane'])
                    M.rootAssembly.rotate(angle=-90.0, 
                                          axisDirection=(0.0, 1.0, 0.0),
                                          axisPoint=(0.0, 0.0, 0.0),
                                          instanceList=('Plane-1', ))
                    M.rootAssembly.translate(vector=(W/2.0, -t, 0),
                                             instanceList=('Plane-1', ))
                    # - Interation property
                    M.ContactProperty('Frictionless')
                    M.interactionProperties['Frictionless'].TangentialBehavior(formulation=FRICTIONLESS)
                    M.interactionProperties['Frictionless'].NormalBehavior(allowSeparation=ON,
                                                                           constraintEnforcementMethod=DEFAULT, 
                                                                           pressureOverclosure=HARD)
                    # - Interaction
                    M.SurfaceToSurfaceContactStd(adjustMethod=NONE,
                                                 clearanceRegion=None,
                                                 createStepName='Initial',
                                                 datumAxis=None, 
                                                 initialClearance=OMIT,
                                                 interactionProperty='Frictionless',
                                                 master=regionToolset.Region(side2Faces=M.rootAssembly.instances['Plane-1'].faces.getSequenceFromMask(mask=('[#f ]', ), )),
                                                 name='Plane-Plate',
                                                 slave=M.rootAssembly.sets['CRACK_FACE_NODES'],
                                                 sliding=FINITE,
                                                 thickness=ON)

                    # Fix boundary conditions
                    M.EncastreBC(createStepName='Initial', 
                                 localCsys=None, name='Fix-RP',
                                 region=regionToolset.Region(referencePoints=(M.rootAssembly.instances['Plane-1'].referencePoints[2],)))

                    # - Delete old ones
                    for i in range(2, 31):
                        del M.steps['Step-{0}'.format(i)]
                    # - Make new one
                    M.StaticStep(name='Step-2',
                                 previous='Step-1',
                                 description='Apply the displacement',
                                 initialInc=0.1,
                                 maxNumInc=1000)
                    allNodes = M.parts['PART-1'].nodes
                    (y, z, delta) = (-t, -S_inner, 1.0e-4)
                    (ymin, ymax, zmin, zmax) = (y-delta, y+delta, z-delta, z+delta)
                    rollerNodes = allNodes.getByBoundingBox(yMin=ymin, yMax=ymax,
                                                            zMin=zmin, zMax=zmax)
                    rollerSet = M.parts['PART-1'].Set(name='Roller',
                                                      nodes=rollerNodes)
                    if enableFullyPlastic:
                        displacement = 2*displacementDict[modelName]
                        elements = M.parts['PART-1'].elements
                        nodes = M.parts['PART-1'].nodes
                        nodes_cf = M.rootAssembly.sets['ContInt-1-Tip'].nodes
                        elements_fp_label = []
                        elements_fp_index = []
                        nodes_fp_label = []
                        nodes_fp_index = []
                        print("Making list of CF node labels")
                        for n in nodes_cf:
                            if n.coordinates[-1] == 0:
                                nodes_fp_label.append(n.label)
                                nodes_fp_index.append(n.label-1)  # for a node index (node label 1 is index 0)
                        nodes_fp_label = list(OrderedDict.fromkeys(nodes_fp_label))  # Remove duplicates, preserve ordering
                        nodes_fp_index = list(OrderedDict.fromkeys(nodes_fp_index))

                        print("Making list of elements along z=0")
                        elements_z0_label = []
                        elements_z0_index = []
                        for e in elements:
                            for n_index in e.connectivity:
                                if (nodes[n_index].coordinates[-1] == 0) and not (e.label in elements_z0_label):
                                    elements_z0_label.append(e.label)
                                    elements_z0_index.append(e.label-1)
                        elements_z0_label = list(OrderedDict.fromkeys(elements_z0_label))  # Remove duplicates, preserve ordering
                        elements_z0_index = list(OrderedDict.fromkeys(elements_z0_index))

                        for ring in range(1, 11):
                            print("Adding elements from ring {0}".format(ring))
                            for n_index in nodes_fp_index:
                                n_label = n_index+1
                                for e_index in elements_z0_index:
                                    if (n_index in elements[e_index].connectivity) and not (elements[e_index].label in elements_fp_label):
                                        elements_fp_label.append(elements[e_index].label)
                                        elements_fp_index.append(elements[e_index].label-1)
                            elements_fp_label = list(OrderedDict.fromkeys(elements_fp_label))  # Remove duplicates, preserve ordering
                            elements_fp_index = list(OrderedDict.fromkeys(elements_fp_index))
                            for e_index in elements_fp_index:
                                for n_index in elements[e_index].connectivity:
                                    if nodes[n_index].coordinates[-1] == 0:
                                        nodes_fp_index.append(n_index)
                                        nodes_fp_label.append(n_index+1)
                                nodes_fp_label = list(OrderedDict.fromkeys(nodes_fp_label))  # Remove duplicates, preserve ordering
                                nodes_fp_index = list(OrderedDict.fromkeys(nodes_fp_index)) 

                        M.rootAssembly.SetFromNodeLabels(name='CHECK_FP_NODES', nodeLabels=(('PART-1-1', nodes_fp_label ), ) )
                        M.rootAssembly.SetFromElementLabels(name='CHECK_FP_ELEMENTS', elementLabels=(('PART-1-1', elements_fp_label ), ) )
                        M.steps['Step-2'].setValues(fullyPlastic='CHECK_FP_ELEMENTS')
                    else:
                        displacement = displacementDict[modelName]
                    M.DisplacementBC(name='Roller',
                                     createStepName='Step-2',
                                     region=M.rootAssembly.instances['PART-1-1'].sets['Roller'],
                                     u2=displacement)

                    # Create job
                    M.rootAssembly.regenerate()
                    mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
                            explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
                            memory=90, memoryUnits=PERCENTAGE, model=modelName,
                            modelPrint=OFF, multiprocessingMode=DEFAULT,
                            name=modelName, nodalOutputPrecision=SINGLE, numCpus=1,
                            numDomains=1, numGPUs=0, queue=None, resultsFormat=ODB, scratch='',
                            type=ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)
