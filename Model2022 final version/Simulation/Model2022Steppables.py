from datetime import datetime
from cc3d.core.PySteppables import *


class Model2022Steppable(SteppableBasePy):
    def __init__(self,frequency=1):
        SteppableBasePy.__init__(self,frequency)
        # References to cell types
        self.typeCollagen = 1
        self.typeTumor = 2
        self.typeNucleus = 3
        self.typeFactin = 4
        
        # Initialize cell references
        self.myTumor = None
        self.myNucleus = None
        
        # Combined target volume of the tumor and factin
        self.totalTargetVolume = 3500
        
        # Other volume atributes of the tumor and nucleus
        self.tumorLambdaVolume = 10
        self.nucleusTargetVolume = 350 # Supposed to be about 10 procent of total target volume
        self.nucleusLambdaVolume = 10
        
        # Percentage of factin 'within' the tumor
        self.newFactinPercentage = 0.3 # Should be between 0.05 and 0.3
        
        # Volume atributes of the newly created factin
        self.newFactinTargetVolume = 15
        self.newFactinLambdaVolume = 10  
        
        # Neighbor order of 'get_pixel_neighbors_based_on_neighbor_order'
        self.BoundaryPixelNeighborOrder = 1
        
        # Clustering counting variable
        self.collagenClusterCounter = 0
        
        # Name & location of the log file
        fileName = 'MyLog'
        filePath = 'C:\\Users\\Femke\\Documents\\Modellen\\eigen modellen\\Model2022 v4\\Simulation' # Adjustable!
        
        # Open log file
        fileLoc = '/'.join([filePath, fileName + '.piff'])
        self.myFile = open(fileLoc, 'w')
        currentTime = datetime.now()
        printTime = currentTime.strftime("%H:%M:%S")
        self.printLog('Steppable initialized on', printTime)

    # Print to log file
    def printLog(self, *args):
        for arg in args:
            self.myFile.write(str(arg) + ' ')
        self.myFile.write('\n')

    # Print array to log file
    def printLogArr(self, arr):
        for elm in arr:
            self.printLog(elm)
    
    # Translates type id to type name
    def typeToName(self, type):
        return ['Medium', 'Collagen', 'Tumor', 'Nucleus', 'Factin'][type]
    

            
    # Create a track with an adjustable narrowing    
    def CollagenCreator(self, upper_x,upper_y,lower_x,lower_y):
        # Fixed dimensions of the start of the track
        fixed_width = 100
        fixed_height = 50
        
        # Create all upper track elements
        collagen_upper_fixed = self.new_cell(self.typeCollagen) # Horizontal element before narrowing
        self.cell_field[self.dim.x-fixed_width:self.dim.x, self.dim.y//2+fixed_height//2, 0] = collagen_upper_fixed
        
        collagen_upper_x = self.new_cell(self.typeCollagen) # Horizontal element of narrowing
        self.cell_field[self.dim.x-fixed_width-upper_x:self.dim.x-fixed_width, self.dim.y//2+fixed_height//2+upper_y, 0] = collagen_upper_x
       
        collagen_upper_x_left = self.new_cell(self.typeCollagen) # Horizontal element after narrowing
        self.cell_field[0:self.dim.x-fixed_width-upper_x, self.dim.y//2+fixed_height//2, 0] = collagen_upper_x_left
        
        # Distinct between a narrowing and a widening
        if upper_y > 0:
            collagen_upper_y = self.new_cell(self.typeCollagen) # Vertical element at the start of the narrowing
            self.cell_field[self.dim.x-fixed_width, self.dim.y//2+fixed_height//2:self.dim.y//2+fixed_height//2+upper_y, 0] = collagen_upper_y
            
            collagen_upper_y_left = self.new_cell(self.typeCollagen) # Vertical element at the end of the narrowing
            self.cell_field[self.dim.x-fixed_width-upper_x, self.dim.y//2+fixed_height//2:self.dim.y//2+fixed_height//2+upper_y, 0] = collagen_upper_y_left
        else:
            collagen_upper_y = self.new_cell(self.typeCollagen) # Vertical element at the start of the narrowing
            self.cell_field[self.dim.x-fixed_width, self.dim.y//2+fixed_height//2+upper_y:self.dim.y//2+fixed_height//2, 0] = collagen_upper_y
            
            collagen_upper_y_left = self.new_cell(self.typeCollagen) # Vertical element at the end of the narrowing
            self.cell_field[self.dim.x-fixed_width-upper_x, self.dim.y//2+fixed_height//2+upper_y:self.dim.y//2+fixed_height//2, 0] = collagen_upper_y_left      
        
        # Assign the same cluster id to all upper track elements
        self.reassign_cluster_id(collagen_upper_fixed, 0)
        self.reassign_cluster_id(collagen_upper_x, 0)
        self.reassign_cluster_id(collagen_upper_x_left, 0)
        self.reassign_cluster_id(collagen_upper_y, 0)
        self.reassign_cluster_id(collagen_upper_y_left, 0)
        
        
        # Create all lower track elements
        collagen_lower_fixed = self.new_cell(self.typeCollagen) # Horizontal element before narrowing
        self.cell_field[self.dim.x-fixed_width:self.dim.x, self.dim.y//2-fixed_height//2, 0] = collagen_lower_fixed
        
        collagen_lower_x = self.new_cell(self.typeCollagen) # Horizontal element of narrowing
        self.cell_field[self.dim.x-fixed_width-lower_x:self.dim.x-fixed_width, self.dim.y//2-fixed_height//2-lower_y, 0] = collagen_lower_x
        
        collagen_lower_x_left = self.new_cell(self.typeCollagen) # Horizontal element after narrowing
        self.cell_field[0:self.dim.x-fixed_width-lower_x, self.dim.y//2-fixed_height//2, 0] = collagen_lower_x_left
        
        # Distinct between a narrowing and a widening
        if lower_y > 0:
            collagen_lower_y = self.new_cell(self.typeCollagen) # Vertical element at the start of the narrowing
            self.cell_field[self.dim.x-fixed_width, self.dim.y//2-fixed_height//2-lower_y:self.dim.y//2-fixed_height//2, 0] = collagen_lower_y
            
            collagen_lower_y_left = self.new_cell(self.typeCollagen) # Vertical element at the end of the narrowing
            self.cell_field[self.dim.x-fixed_width-lower_x, self.dim.y//2-fixed_height//2-lower_y:self.dim.y//2-fixed_height//2, 0] = collagen_lower_y_left
        else:
            collagen_lower_y = self.new_cell(self.typeCollagen) # Vertical element at the start of the narrowing
            self.cell_field[self.dim.x-fixed_width, self.dim.y//2-fixed_height//2:self.dim.y//2-fixed_height//2-lower_y, 0] = collagen_lower_y
            
            collagen_lower_y_left = self.new_cell(self.typeCollagen) # Vertical element at the end of the narrowing
            self.cell_field[self.dim.x-fixed_width-lower_x, self.dim.y//2-fixed_height//2:self.dim.y//2-fixed_height//2-lower_y, 0] = collagen_lower_y_left      
        
        # Assign the same cluster id to all lower track elements
        self.reassign_cluster_id(collagen_lower_fixed, 1)
        self.reassign_cluster_id(collagen_lower_x, 1)
        self.reassign_cluster_id(collagen_lower_x_left, 1)
        self.reassign_cluster_id(collagen_lower_y, 1)
        self.reassign_cluster_id(collagen_lower_y_left, 1)
        
    # Create a track that is diagonal at the end of the narrowing    
    def DiagonalCollagenCreator(self):
        # Set all dimensions
        fixed_width = 100
        fixed_height = 50
        upper_x = 100
        upper_y = -10
        lower_x = 100
        lower_y = -10
        
        # Create all upper track elements
        collagen_upper_fixed = self.new_cell(self.typeCollagen) # Horizontal element before narrowing
        self.cell_field[self.dim.x-fixed_width:self.dim.x, self.dim.y//2+fixed_height//2, 0] = collagen_upper_fixed
        
        collagen_upper_x = self.new_cell(self.typeCollagen) # Horizontal element of narrowing
        self.cell_field[self.dim.x-fixed_width-upper_x:self.dim.x-fixed_width, self.dim.y//2+fixed_height//2+upper_y, 0] = collagen_upper_x
       
        collagen_upper_x_left = self.new_cell(self.typeCollagen) # Horizontal element after narrowing
        self.cell_field[0:self.dim.x-fixed_width-upper_x-10, self.dim.y//2+fixed_height//2, 0] = collagen_upper_x_left
        
        
        collagen_upper_y = self.new_cell(self.typeCollagen) # Vertical element at the start of the narrowing
        self.cell_field[self.dim.x-fixed_width, self.dim.y//2+fixed_height//2+upper_y:self.dim.y//2+fixed_height//2, 0] = collagen_upper_y
        
        # Create all steps at the end of the narrowing
        collagen_upper_step1 = self.new_cell(self.typeCollagen) 
        self.cell_field[self.dim.x-fixed_width-upper_x-2:self.dim.x-fixed_width-upper_x,self.dim.y//2+fixed_height//2+upper_y:self.dim.y//2+fixed_height//2+upper_y+2,0] = collagen_upper_step1
        
        collagen_upper_step2 = self.new_cell(self.typeCollagen)
        self.cell_field[self.dim.x-fixed_width-upper_x-4:self.dim.x-fixed_width-upper_x-2,self.dim.y//2+fixed_height//2+upper_y+2:self.dim.y//2+fixed_height//2+upper_y+4,0] = collagen_upper_step2
        
        collagen_upper_step3 = self.new_cell(self.typeCollagen)
        self.cell_field[self.dim.x-fixed_width-upper_x-6:self.dim.x-fixed_width-upper_x-4,self.dim.y//2+fixed_height//2+upper_y+4:self.dim.y//2+fixed_height//2+upper_y+6,0] = collagen_upper_step3
        
        collagen_upper_step4 = self.new_cell(self.typeCollagen)
        self.cell_field[self.dim.x-fixed_width-upper_x-8:self.dim.x-fixed_width-upper_x-6,self.dim.y//2+fixed_height//2+upper_y+6:self.dim.y//2+fixed_height//2+upper_y+8,0] = collagen_upper_step4
        
        collagen_upper_step5 = self.new_cell(self.typeCollagen)
        self.cell_field[self.dim.x-fixed_width-upper_x-10:self.dim.x-fixed_width-upper_x-6,self.dim.y//2+fixed_height//2+upper_y+8:self.dim.y//2+fixed_height//2+upper_y+10,0] = collagen_upper_step5
        
        # Assign the same cluster id to all upper track elements
        self.reassign_cluster_id(collagen_upper_fixed, 0)
        self.reassign_cluster_id(collagen_upper_x, 0)
        self.reassign_cluster_id(collagen_upper_x_left, 0)
        self.reassign_cluster_id(collagen_upper_y, 0)
        self.reassign_cluster_id(collagen_upper_step1, 0)
        self.reassign_cluster_id(collagen_upper_step2, 0)
        self.reassign_cluster_id(collagen_upper_step3, 0)
        self.reassign_cluster_id(collagen_upper_step4, 0)
        self.reassign_cluster_id(collagen_upper_step5, 0)
        
        
        # Create all lower track elements
        collagen_lower_fixed = self.new_cell(self.typeCollagen) # Horizontal element before narrowing
        self.cell_field[self.dim.x-fixed_width:self.dim.x, self.dim.y//2-fixed_height//2, 0] = collagen_lower_fixed
        
        collagen_lower_x = self.new_cell(self.typeCollagen) # Horizontal element of narrowing
        self.cell_field[self.dim.x-fixed_width-lower_x:self.dim.x-fixed_width, self.dim.y//2-fixed_height//2-lower_y, 0] = collagen_lower_x
        
        collagen_lower_x_left = self.new_cell(self.typeCollagen) # Horizontal element after narrowing
        self.cell_field[0:self.dim.x-fixed_width-lower_x-10, self.dim.y//2-fixed_height//2, 0] = collagen_lower_x_left
        
        
        collagen_lower_y = self.new_cell(self.typeCollagen) # Vertical element at the start of the narrowing
        self.cell_field[self.dim.x-fixed_width, self.dim.y//2-fixed_height//2:self.dim.y//2-fixed_height//2-lower_y, 0] = collagen_lower_y
        
        # Create all steps at the end of the narrowing 
        collagen_lower_step1 = self.new_cell(self.typeCollagen)
        self.cell_field[self.dim.x-fixed_width-lower_x-2:self.dim.x-fixed_width-lower_x,self.dim.y//2-fixed_height//2-lower_y-2:self.dim.y//2-fixed_height//2-lower_y,0] = collagen_lower_step1
        
        collagen_lower_step2 = self.new_cell(self.typeCollagen)
        self.cell_field[self.dim.x-fixed_width-lower_x-4:self.dim.x-fixed_width-lower_x-2,self.dim.y//2-fixed_height//2-lower_y-4:self.dim.y//2-fixed_height//2-lower_y-2,0] = collagen_lower_step2
        
        collagen_lower_step3 = self.new_cell(self.typeCollagen)
        self.cell_field[self.dim.x-fixed_width-lower_x-6:self.dim.x-fixed_width-lower_x-4,self.dim.y//2-fixed_height//2-lower_y-6:self.dim.y//2-fixed_height//2-lower_y-4,0] = collagen_lower_step3
        
        collagen_lower_step4 = self.new_cell(self.typeCollagen)
        self.cell_field[self.dim.x-fixed_width-lower_x-8:self.dim.x-fixed_width-lower_x-6,self.dim.y//2-fixed_height//2-lower_y-8:self.dim.y//2-fixed_height//2-lower_y-6,0] = collagen_lower_step4
        
        collagen_lower_step5 = self.new_cell(self.typeCollagen)
        self.cell_field[self.dim.x-fixed_width-lower_x-10:self.dim.x-fixed_width-lower_x-6,self.dim.y//2-fixed_height//2-lower_y-10:self.dim.y//2-fixed_height//2-lower_y-8,0] = collagen_lower_step5
         
        # Assign the same cluster id to all lower track elements
        self.reassign_cluster_id(collagen_lower_fixed, 1)
        self.reassign_cluster_id(collagen_lower_x, 1)
        self.reassign_cluster_id(collagen_lower_x_left, 1)
        self.reassign_cluster_id(collagen_lower_y, 1)
        self.reassign_cluster_id(collagen_lower_step1, 1)
        self.reassign_cluster_id(collagen_lower_step2, 1)
        self.reassign_cluster_id(collagen_lower_step3, 1)
        self.reassign_cluster_id(collagen_lower_step4, 1)
        self.reassign_cluster_id(collagen_lower_step5, 1)
        
            
    # Create various tracks with different dimensions
    def CollagenTrack(self, tracknumber):
        # y values may not exceed 55
        # x values may not exceed 300
        # Input are the dimensions of the narrowing= upper_x,upper_y,lower_x,lower_y 
        
        # Straight tube
        if tracknumber == 1: 
            self.CollagenCreator(0,0,0,0)
        # Upper indent
        elif tracknumber == 2:
            self.CollagenCreator(100,-30,0,0)
        # Lower indent
        elif tracknumber == 3:
            self.CollagenCreator(0,0,100,-12)
        # Upper and lower indent
        elif tracknumber == 4:
            self.CollagenCreator(150,-10,150,-10)
        # Upper and lower indent shifted
        elif tracknumber == 5:
            self.CollagenCreator(100,-12,150,-12)
        # Sharp upper and lower indent
        elif tracknumber == 6:
            self.CollagenCreator(10,-15,10,-15)
        # Broadening
        elif tracknumber == 7:
            self.CollagenCreator(300,5,300,5)
        # Diagonal version of track 4
        elif tracknumber == 8:
            self.DiagonalCollagenCreator()
        
    
    def start(self):     
        
        self.CollagenTrack(1) # Adjustable!
            
        # Note the tumor and nucleus need to be created after the collagen,
        # to prevent the cell cluster and collagen cluster to coincide
        
        # Create the tumor
        tumor = self.new_cell(self.typeTumor)
        tumor.targetVolume = self.totalTargetVolume
        tumor.lambdaVolume = self.tumorLambdaVolume
        self.cell_field[self.dim.x-50, self.dim.y//2, 0] = tumor
        self.myTumor = tumor
        
        # Create the nucleus
        nucleus = self.new_cell(self.typeNucleus)
        nucleus.targetVolume = self.nucleusTargetVolume
        nucleus.lambdaVolume = self.nucleusLambdaVolume
        self.cell_field[self.dim.x-40, self.dim.y//2, 0] = nucleus
        self.myNucleus = nucleus
        
        # An example of how a graph would be added
        #self.plot_win = self.add_new_plot_window(
        #title='Relative distance of the Center of Mass of the Tumor and Nucleus',
        #x_axis_title='MonteCarlo Step*10 (MCS)',
        #y_axis_title='Relative distance',
        #x_scale_type='linear',
        #y_scale_type='linear',
        #grid=True  )
        #self.plot_win.add_plot("Relative distance", style='Dots', color='red', size=5)


    # Compute the total target volume of a given cell type
    def totalTargetVolumeByCell(self, type):
        totalVolume = 0
        for cell in self.cell_list_by_type(type):
            totalVolume += cell.targetVolume
        return totalVolume    

    # Compute the total target volume of a given cell type within a cluster
    def totalTargetVolumeByCellGivenClusters(self, type, clusterSet):
        totalVolumeArray = [0] * (len(clusterSet))
        for cell in self.cell_list_by_type(type):
            totalVolumeArray[cell.clusterId] += cell.targetVolume
        return totalVolumeArray

    # Swap collagen clusters
    def swapCollagenClusters(self, oldCluster, newCluster):
        for collagen in self.cell_list_by_type(self.typeCollagen):
            if collagen.clusterId == oldCluster:
                self.reassign_cluster_id(collagen, newCluster)
            elif collagen.clusterId == newCluster:
                self.reassign_cluster_id(collagen, oldCluster)

    def step(self,mcs):
        # Check if my tumor is created
        if self.myTumor:
            # Determine the total target volumes
            totalTargetVolumeFactin = self.totalTargetVolumeByCell(self.typeFactin)
            totalTargetVolumeMyCell = self.myTumor.targetVolume + totalTargetVolumeFactin
            
            # Determine the (in)direct collagen neighbors     
            collagenSet = set()
            for neighbor, commonSurfaceArea in self.get_cell_neighbor_data_list(self.myTumor):
                # Check if neigbor exists
                if neighbor:
                    # Add the neighbor to the collagen set
                    if neighbor.type == self.typeCollagen:
                        collagenSet.add(neighbor.clusterId)
                    if neighbor.type == self.typeFactin:
                        for neighborOrder2, commonSurfaceAreaOrder2 in self.get_cell_neighbor_data_list(neighbor):
                            # Check if neigbor of order 2 exists
                            if neighborOrder2:
                                # Add the neighbor of order 2 to the collagen set
                                if neighborOrder2.type == self.typeCollagen:
                                    collagenSet.add(neighborOrder2.clusterId)
            
            # Complete the set with previous collagen neighbors
            for factin in self.cell_list_by_type(self.typeFactin):
                collagenSet.add(factin.clusterId)
            
            # Determine the total factin target volumes within the clusters
            totalFactinInCluster = self.totalTargetVolumeByCellGivenClusters(self.typeFactin, collagenSet)
            
            # Correct cluster numbering
            singeltonZero = set()
            singeltonZero.add(0)
            singeltonOne = set()
            singeltonOne.add(1)
            if collagenSet == singeltonOne:
                # Swap collagen clusters
                self.swapCollagenClusters(0, 1)
                # Reconstruct the collagen set after swapping
                collagenSet = singeltonZero
                # Debug statement
                self.printLog('Collagen cluster order is swapped')
            
            # Iterate over the collagen clusters
            numberOfCollagenClusters = len(collagenSet)
            for collagenCluster in collagenSet:
                # Determine if there is still intercelluar factin left
                allowedFactinInCluster = self.totalTargetVolume * (self.newFactinPercentage / numberOfCollagenClusters)
                excessFactin = totalFactinInCluster[collagenCluster] - allowedFactinInCluster
                
                # Add new Factin
                if excessFactin < 0:
                    # Determine the boundary points
                    for point in self.get_cell_boundary_pixel_list(self.myTumor):
                        for neighborPoint in self.get_pixel_neighbors_based_on_neighbor_order(point.pixel, self.BoundaryPixelNeighborOrder):
                            # Get the neighbor cell
                            neighbor = self.cell_field[neighborPoint.pt.x, neighborPoint.pt.y, 0]
                            # Check if neighbor exists
                            if neighbor:
                                # Check if the neighbor is collagen and within the current cluster- 
                                if neighbor.type == self.typeCollagen and neighbor.clusterId == collagenCluster:
                                    # Create factin on the current location
                                    newFactin = self.new_cell(self.typeFactin)
                                    newFactin.targetVolume = self.newFactinTargetVolume
                                    newFactin.lambdaVolume = self.newFactinLambdaVolume
                                    self.reassign_cluster_id(newFactin, collagenCluster)
                                    self.cell_field[int(point.pixel.x), int(point.pixel.y),0] = newFactin
                # Reduce the amount of Factin in this cluster
                else:
                    for factin in self.cell_list_by_type(self.typeFactin):
                        # Check if the factin is within the cluster
                        if factin.clusterId == collagenCluster:
                            # Determine the removable target volume
                            removableTargetVolume = min(factin.targetVolume, excessFactin)
                            excessFactin -= removableTargetVolume
                            # Reduce the target volume if possible
                            if excessFactin < 0:
                                break
                            else:
                                factin.targetVolume -= removableTargetVolume   
            
            # Update the total target volumes
            totalTargetVolumeFactin = self.totalTargetVolumeByCell(self.typeFactin)
            totalTargetVolumeMyCell = self.myTumor.targetVolume + totalTargetVolumeFactin
            total = totalTargetVolumeFactin + totalTargetVolumeMyCell
            
            # Manage the total target volume
            self.myTumor.targetVolume = self.totalTargetVolume - totalTargetVolumeFactin
            
            # Merge the newly created factin
            for factin in self.cell_list_by_type(self.typeFactin):
                for factinNeighbor, commonSurfaceArea in self.get_cell_neighbor_data_list(factin):
                    # Check if neighbor exists
                    if factinNeighbor:
                        # Check if neighbor is factin and not itself
                        if factinNeighbor.type == self.typeFactin and factin.id != factinNeighbor.id:
                            # Check if neighbor and factin are within the same cluster
                            if factinNeighbor.clusterId == factin.clusterId:
                                # Check if the target volume is non-zero
                                if factin.targetVolume != 0:
                                    factin.targetVolume += factinNeighbor.targetVolume
                                    factinNeighbor.targetVolume = 0
            
            # Delete factin with target volume 0
            for factin in self.cell_list_by_type(self.typeFactin):
                if factin.targetVolume == 0:
                    self.delete_cell(factin)   
        
        # An example of how a graph would be created (continued)
        # Check if my nucleus is created
        #if self.myNucleus:
        #    if mcs % 10 == 0:
        #        self.plot_win.add_data_point('Relative distance', mcs/10, self.myTumor.xCOM-self.myNucleus.xCOM)


    def on_stop(self):
        # Debug statement
        self.printLog()
        self.printLog('Cell list when stopped')
        self.printLog('type, cluster id, target volume')
        for cell in self.cell_list:
            self.printLog(self.typeToName(cell.type), cell.clusterId, cell.targetVolume)
        self.printLog()
        
        self.printLog('Steppable is stopped')
        self.myFile.close()
        return

    def finish(self):
        self.printLog('Steppable is finished')
        self.myFile.close()
        return
