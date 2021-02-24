import numpy as np 
import operator as op

one = 1 
zero = 0
 
class GlobalAlignment:
    def __init__(self , sequence1 , sequence2 ): # made a constructor encompassing the 2 input sequences for our case
        self.sequence1 = sequence1; 
        self.sequence2 = sequence2; 
        
    def scoring(self):
        print(' ') 
        
    #----SCORING CONDITIONS----# 
    print("Input the SCORING CONDITIONS:" )
    print(" ")
    match = int(input("Input the MATCH :" )) #taking the respective inputs for all the three scoring condtions match 
    mismatch = int(input("Input the MISMATCH : ")) # input for mismatch 
    gap_Penalty = int(input("Input the GAP PENALTY :" )) # input for gap penalty 
    # match = 2
    # mismatch = -1
    # gap_Penalty = -1
    def globalAlignmentScore(self): 
        sequence_1 = list(self.sequence1) #converts the given string sequence into a list in a format of characters split with respect to commas 
        sequence_2 = list(self.sequence2)
        rows = sequence_1[:] # ['A', 'T', 'C', 'A', 'G', 'A', 'G', 'T', 'A']
        cols = sequence_2[:] # ['T', 'T', 'C', 'A', 'G', 'T', 'A']

        rows.insert(0,0) #this inserts gap element at the 0th index of the 1st row having no value 
        cols.insert(0,0) #this inserts gap element at the 0th index of the 1st column having no value 
        # we do this for the purpose of having a gap element at the beginning of the Matrix
        
        # ----- DP MATRIX CREATION ----- #
        DP_Matrix = [] #initialising the DP MATRIX 
        DP_Matrix = np.zeros((np.size(rows),np.size(cols))) #we are making a DP_Matrix having only 0s in all the postions
        count_Column = 0 #we assign an intger count_Column to keep track of the previous iteration's value with respect to the subtraction made with the Mismatch penalty 
        x = 0 
        while x in range(np.size(cols)): #looping over the range of the size of column 
            DP_Matrix[0][x] = count_Column #assigning values to the 1st column 
            count_Column = op.add(count_Column,GlobalAlignment.mismatch) #performs the addition between the mismatch penalty and the count_Column
            x = x + 1 # loops through all indexes until the length of the 1st column is over 
        count_Row = 0 #we assign an intger count_Row to keep track of the previous iteration's value with respect to the subtraction made with the Mismatch penalty 
        x = 0 
        while x in range(np.size(rows)): 
            DP_Matrix[x][0] = count_Row #assigning values to the 1st row 
            count_Row = op.add(count_Row,GlobalAlignment.gap_Penalty) #here we are assigning elements of the first row values with respect to the mismatch penalty globalScore of the previous index 
            x = x + 1  # loops through all indexes until the length of the 1st row is over 
            
        # the DP_Matrix[0][0] = 0 because we take into account the count_Row = 0 and count_Column = 0

        row = 1 
        while row in range(1, np.size(rows)): #loops through all the iterations from 2nd row onwards in this loop as we have already assigned values to the 1st column and 1st row
            col = 1 
            while col in range(1, np.size(cols)):
                temp_List = [] # we create a temporary List to get the max of all the values we compute for a particular i and j in the DP_Matrix
                #this condition checks if the previous characters in 2 DNA Sequences are equal or not.
                if op.ne(sequence_1[row-1],sequence_2[col-1]):  # if unequal, the value at D[row][column] becomes the addition of DP_Matrix[row-1][col-1] and the match globalScore 
                    temp_List.extend([DP_Matrix[row-1][col-1] + GlobalAlignment.mismatch]) 
                else: # if equal, the value at D[row][column] becomes the addition of DP_Matrix[row-1][col-1] and the match globalScore 
                    temp_List.extend([DP_Matrix[row-1][col-1] + GlobalAlignment.match])
                temp_List.append(DP_Matrix[row-1][col] + GlobalAlignment.gap_Penalty) #by default we compute the addition of the value to the left of the DP_matrix and gap penalty  
                temp_List.append(DP_Matrix[row][col-1] + GlobalAlignment.gap_Penalty) #by default we compute the addition of the value to the above of the DP_matrix and gap penalty   
                DP_Matrix[row][col] = int(max(temp_List)) # we compute the maximum of all the values we obtain in the temporary List
                col = col + 1
            row = row + 1 
            
        p2 = PrintGrid(sequence_1, sequence_2 ,DP_Matrix); 
        p2.printDPMatrix();
        
class PrintGrid: #this is my PrintGrid Class
    def __init__( self , sequence_1 , sequence_2 , DP_Matrix ): #constructor calling the 2 sequences and the DP matrix
        self.sequence_1 = sequence_1; 
        self.sequence_2 = sequence_2;
        self.DP_Matrix = DP_Matrix;
    
    def printDPMatrix(self):
        
        
        max = 0; #this is the max globalScore of this Global Alignment Matrix, which is assigned to 0 for the sake of keeping track
        print('DP Matrix')
        for x in self.DP_Matrix: #loops through all the rows of the given DP_Matrix
            print(x) #prints each and every row of the DP Matrix 
            for y in x : #loops through each and every element of the row x 
                if op.gt(y ,max): #checks if the element y from the row is greater than the max globalScore 
                    max = y #if the element y is greater, value y is assigned to Max 
                    continue 
                else:
                    continue
        
        print("Max globalScore = ", int(max)) #prints the max element of the given sequence
        print(" ") 
        p3 = BackTracking(self.sequence_1 , self.sequence_2 , self.DP_Matrix , max )
        p3.backtracking();

final_1st = [] # assigned an empty list where we add all the possible aligned sequences of the 1st Sequence 
final_2nd = [] # assigned an empty list where we add all the possible aligned sequences of the 2nd Sequence  
astring = '' 
bstring = ''

class BackTracking: #
    def __init__(self , sequence_1 , sequence_2 , DP_Matrix , max): #constructor containing the 2 seqeunces , DP matrix and the maximum element
        self.sequence_1 = sequence_1;
        self.sequence_2 = sequence_2;
        self.DP_Matrix = DP_Matrix ; 
        self.max = max ;  
        
    def backtracking(self):
        i = np.size(self.sequence_1) # here we take the final ith index of the last element of the dp matrix , taken with respect to the rows 
        j = np.size(self.sequence_2) # here we take the final jth index of the last element of the dp matrix , taken with respect to the columns
        k = max(i,j)
        stringSet = [('', '')] 
        maxIndexSet = [(i, j)] 
         
        def tempList(i, j , k): 
            # initiated a temporary list where we keep track of the several possibilities the arrow can be pointed towards once we begin back tracking 
            tmpCheckList = []
            if op.eq(op.add(self.DP_Matrix[i-one][j-one],GlobalAlignment.match),(self.DP_Matrix[i][j])):
                if op.eq(self.sequence_1[i-one],self.sequence_2[j-one]):
                    tmpCheckList.extend([1]) #we add 1 to this list if the ma#list of possible final aligned forms of sequence 1tch condition is proved to be correct       
            if op.eq(op.add(self.DP_Matrix[i-1][j-one],GlobalAlignment.mismatch),self.DP_Matrix[i][j]):
                if(one>zero):
                    tmpCheckList.extend([one]) #we add 1 to this list if the mismatch condition is proved to be correct 
            if op.eq(op.add(self.DP_Matrix[i-one][j*one],GlobalAlignment.gap_Penalty),self.DP_Matrix[i][j]):
                if one>zero: 
                    tmpCheckList.extend([2]) #we add 2 to this list if the gap pentalty condition is proved to be correct with respect to the left side 
            if op.eq(op.add(self.DP_Matrix[i][j-one],GlobalAlignment.gap_Penalty),self.DP_Matrix[i][j]): 
                tmpCheckList.extend([3])  #we add 3 to this list if the gap pentalty condition is proved to be correct with respect to the top side     
            return tmpCheckList
        def globalScore(i, j, tmpCheckList, sequenceA, sequenceB , maxIndexSet): 
            a = i ;
            b = j ; 
            while op.ne(np.size(stringSet),0): 
                if(np.size(stringSet)>0): 
                    if op.gt(np.size(tmpCheckList),1): # more than 1 possibilities for having an alignment
                        maxIndexSet[0:0] = [(i, j )] 
                        temporary_List[0:0] = [tmpCheckList[1:]]
                        stringSet[0:0] = [(sequenceA,sequenceB)]
                    if op.gt(np.size(tmpCheckList),0): 
                        # we call this conditional to check if there are more than 1 possibility of having an alignment
                        if op.eq(tmpCheckList[0*1],1): 
                            sequenceA = op.add(self.sequence_1[a-1],sequenceA) 
                            sequenceB = op.add(self.sequence_2[b-1], sequenceB)
                            a = op.sub(a,1)
                            b = op.sub(b,1)
                            tmpCheckList.pop(0) #here we get rid of the current possibility as we already accounted for it  
                        elif op.eq(tmpCheckList[0*1],2):
                            sequenceA = self.sequence_1[a-1] + sequenceA
                            sequenceB = "_" + sequenceB 
                            a = op.sub(a,1)
                            tmpCheckList.pop(0) #here we get rid of the current possibility as we already accounted for it
                        elif op.eq(tmpCheckList[0*1],3):
                            sequenceA = "_" + sequenceA 
                            sequenceB = self.sequence_2[b-1] + sequenceB
                            b = op.sub(b,1) 
                            tmpCheckList.pop(0) #here we get rid of the current possibility as we already accounted for it
                        globalScore(a, b, tempList(a,b,0), sequenceA, sequenceB,maxIndexSet) 
                        #we recursively call this function again and again  until we end up with no more possibilities
                    else: 
                        if op.eq(a,0): # this conditional means we reached the TopLeft DP_Matrix[0][0] 
                            if op.eq(b,0) :
                                final_1st.append(sequenceA) # we finally add the given sequences to the final list for output later on 
                                final_2nd.append(sequenceB)
                        if op.ne(np.size(stringSet), 0): #keeping in mind more than 1 possibility, we recursively keep calling until we are done with everything
                            tempstringa = stringSet[0][0]
                            tempstringb = stringSet[0][1]
                            tempi = maxIndexSet[0][0]
                            tempj = maxIndexSet[0][1]
                            tempcheckarr = temporary_List[0]
                            deleteOperations(stringSet , maxIndexSet , temporary_List)
                            globalScore (tempi, tempj, tempcheckarr, tempstringa, tempstringb , maxIndexSet) 
                else: 
                    break            
        
        temporary_List = [tempList(i, j ,k)]
                
        def deleteOperations(stringSet , maxIndexSet, temporary_List):
            del (stringSet[0] , maxIndexSet[0] , temporary_List[0]) 
        globalScore(i, j, temporary_List[0], astring, bstring , maxIndexSet)
        count = 1;
        x = 0
        while x in range(np.size(final_1st)): 
            print ("Count =" , count  )#Number of possibilities 
            print(final_1st[x]) #prints the 1st 
            print(final_2nd[x]) #prints the 2nd Sequence Alignment 
            print ('Score =' , int(self.max))  #prints the globalScore 
            print(" ") 
            x = op.add(x,1)
            count = op.add(count,1)
p1 = GlobalAlignment('ATCAGAGTA', 'TTCAGTA')
p1.scoring()
p1.globalAlignmentScore()