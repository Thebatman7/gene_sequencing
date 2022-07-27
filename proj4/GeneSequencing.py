#!/usr/bin/python3

from ctypes import alignment
from socket import if_indextoname
from sunau import AUDIO_FILE_ENCODING_LINEAR_16
from which_pyqt import PYQT_VER
if PYQT_VER == 'PYQT5':
	from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
	from PyQt4.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT6':
	from PyQt6.QtCore import QLineF, QPointF
else:
	raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))

import math
import time
import random

# Used to compute the bandwidth for banded version
MAXINDELS = 3

# Used to implement Needleman-Wunsch scoring
MATCH = -3
INDEL = 5
SUB = 1

class GeneSequencing:

	def __init__( self ):
		pass

# This is the method called by the GUI.  _seq1_ and _seq2_ are two sequences to be aligned, _banded_ is a boolean that tells
# you whether you should compute a banded alignment or full alignment, and _align_length_ tells you
# how many base pairs to use in computing the alignment
#                     i     j              
	def align(self, seq1, seq2, banded, align_length):
		self.banded = banded
		self.MaxCharactersToAlign = align_length
		# print("B Length of seq1 %s  seque1: %s" % (len(seq1), seq1))
		# print("B Length of seq2 %s  seque2: %s" % (len(seq2), seq2))
		seq1 = seq1[:align_length]
		seq2 = seq2[:align_length]
		# print("Align Length : %s" % self.MaxCharactersToAlign)
		# print("B Length of seq1 %s  seque1: %s" % (len(seq1), seq1))
		# print("B Length of seq2 %s  seque2: %s" % (len(seq2), seq2))

		
		if(banded):
			result = self.bandedAlignment(seq1, seq2)
		else:
			result = self.unbandedAlignment(seq1, seq2)
		
		score = result[0]
		alignment1 = result[1]
		alignment2 = result[2]

###################################################################################################
# your code should replace these three statements and populate the three variables: score, alignment1 and alignment2
		# score = random.random()*100;
		# alignment1 = 'abc-easy  DEBUG:({} chars,align_len={}{})'.format(
		# 	len(seq1), align_length, ',BANDED' if banded else '')
		# alignment2 = 'as-123--  DEBUG:({} chars,align_len={}{})'.format(
		# 	len(seq2), align_length, ',BANDED' if banded else '')
###################################################################################################

		return {'align_cost':score, 'seqi_first100':alignment1, 'seqj_first100':alignment2}

	#  E(i, j) = min[csub or match + E(i-1, j-1),  cinsert + E(i, j-1),  cdelete + E(i-1, j)]       
	#  diff():  csub = 1   or  cmatch = -3;  cindel = cinsert = cdelete = 5
	def unbandedAlignment(self, seq1, seq2):
		#NeedleMan-Wunch costs:
		c_insert = 5
		c_delete = 5
		c_sub = 1
		c_match = -3
		#self.printInfo(seq1, seq2, align_length)
		matrix = [[]]
		#print(matrix)
		rows, cols = (len(seq1) + 1, len(seq2) + 1)
		#print("rows: %s cols %s" % (rows, cols))
		matrix = [[0 for i in range(cols)] for j in range(rows)]
		#self.printMatrix(matrix)

		pointers = [[]]
		pointers = [[(0, 0) for i in range(cols)] for j in range(rows)]
		for i in range(1): #O(length i)
			for j in range(1, cols):
				pointers[i][j] = (i, j-1)
		for j in range(1): #O(length j)
			for i in range(1, rows):
				pointers[i][j] = (i-1, j)
		#print("Printing pointers matrix:")
		#self.printPointerMatrix(pointers)	

		
		for i in range(1): #O(length i)
			#for j in range(1, cols):
			for j in range(cols):
				matrix[i][j] = j * c_insert
				#matrix[i][j] = seq2[j - 1]
		for j in range(1): #O(length j)
			#for i in range(1, rows):
			for i in range(rows):
				matrix[i][j] = i * c_delete	
				#matrix[i][j] = seq1[i - 1]
		#self.printMatrix(matrix)
		#print("we have printed table")

		diagonal = 0
		left = 0
		top = 0
		options = [] #left, top, diagonal 
		for i in range(1, rows): 
			for j in range(1, cols):
				if(seq1[i-1] == seq2[j-1]):
					diagonal = matrix[i-1][j-1] - 3
				else:
					diagonal = matrix[i-1][j-1] + 1
				left = matrix[i][j-1] + 5
				top = matrix[i-1][j] + 5

				#break ties with the following preference order: left, top, diagonal.
				options.append(left)
				options.append(top)
				options.append(diagonal)
				min_v = min(options)#O(3) = O(1)

				matrix[i][j] = min_v
				index = options.index(min_v)#O(3) = O(1)
				
				if(index == 0): #left
					pointers[i][j] = (i,j-1)
				if(index == 1): #top
					pointers[i][j] = (i-1,j)
				if(index == 2): #diagonal
					pointers[i][j] = (i-1,j-1)

				options = []

		#self.printMatrix(matrix)
		#self.printPointerMatrix(pointers)
		#print("")

		score = matrix[rows-1][cols-1]
		origin = (0,0)
		tuple = pointers[rows - 1][cols - 1]
		dist = (rows-1, cols-1) #goal
		align1 = ""
		align2 = ""
		# print(seq1)
		# print(len(seq1))
		# print(seq2)
		# print(len(seq2))
		# print(tuple)
		# print("This is the value of i %s" % tuple[0])
		# print("This is the value of j %s" % tuple[1])
		# print("This is the char value of sequ1 at index %s, %s" % (tuple[0], seq1[dist[0]-1]))
		# print("This is the char value of sequ2 at index %s, %s" % (tuple[1], seq2[dist[1]-1]))
		while(origin != dist):
			i_index = dist[0] - tuple[0]
			j_index = dist[1] - tuple[1]

			#how are we supposed to modify the the strings ??? Score
			if(i_index == 1 and j_index == 1): # diagonal 
				align1 = seq1[dist[0] - 1] + align1
				align2 = seq2[dist[1] - 1] + align2
			
			if(i_index == 1 and j_index == 0): # top
				align1 = seq1[dist[0] - 1] + align1
				align2 = "-" + align2
				
		
			if(i_index == 0 and j_index == 1): # left
				align1 = "-" + align1
				align2 = seq2[dist[1] - 1] + align2

			dist = tuple
			tuple = pointers[dist[0]][dist[1]]

		align1 = align1[:100]
		align2 = align2[:100]
		# print(align1)
		# print(align2)
		return score, align1, align2



	#Second implementation
	def bandedAlignment(self, seq1, seq2):
		#NeedleMan-Wunch costs:
		c_insert = 5
		c_delete = 5
		c_sub = 1
		c_match = -3
		bandWith = (2 * 3) + 1 # 2*d+1
		if((len(seq2) - len(seq1)) > 3):
			return float('inf'), "No Alignment Possible", "No Alignment Possible"
		else :
			rows = len(seq1) + 1
			cols = bandWith
			#self.printBandedInf(seq1, seq2, rows, cols)
			table = [[float("inf") for i in range(cols)] for j in range(rows)] # inf middle first row

			pointers = [[(float("inf"), float("inf")) for i in range(cols)] for j in range(rows)]

			#self.printTable(table)
			#self.printT(table)
			#range(startindex, stopindex (not inclusive), step)
			for j in range(3, -1, -1): # O(4)  # /
				table[3-j][j] = (3-j)*5 
				if(j != 0):
					pointers[j][3-j] = (j-1, (3-j)+1)

			for i in range(1): # O(1) # __
				for j in range(4): # O(4)
					table[i][j+3] = j * 5
					if(j != 0):
						pointers[i][j+3] = (0, j+2)
			pointers[0][3] = (0,0)
			#self.printPointerMatrix(pointers)


			start = 4
			table_min_val = float("inf")
			t_min_val_location = (0,0)
			
			diagonal = float("inf")
			left = float("inf")
			top = float("inf")
			options = [] # left, top, diagonal 
			for i in range(1, rows): # O(n-1) = O(n)
				if(start > 0):
					start -= 1
				for j in range(start, cols, 1): # O(7) = O(1)
					index_of_seq2 = (i + j) - 4
					if(index_of_seq2 >= len(seq2)):
						break

					if(seq1[i-1] == seq2[index_of_seq2]):
						#print("char in seq1 " + seq1[i-1] + " " "char in seq2 " + seq2[index_of_seq2])
						diagonal = table[i-1][j] - 3
					else:
						diagonal = table[i-1][j] + 1
					
					if(j != cols - 1): # right limit
						top = table[i-1][j+1] + 5

					if(j != 0): # left limit
						left = table[i][j-1] + 5
					
					options.append(left)
					options.append(top)
					options.append(diagonal)
					min_val = min(options) # O(3) = O(1)

					table[i][j] = min_val
			

					if(min_val < table_min_val): # O(1)
						table_min_val = min_val # O(1)
						t_min_val_location = (i, j)
					

					# tracking back pointer
					index = options.index(min_val) # O(3) = O(1)
					
					if(index == 0): #left
						pointers[i][j] = (i,j-1)
					if(index == 1): #top
						pointers[i][j] = (i-1,j+1)
					if(index == 2): #diagonal
						pointers[i][j] = (i-1,j)

					options = []
					diagonal = float("inf")
					left = float("inf")
					top = float("inf")


			j_ind = (len(seq2) + 4) - rows
			table_min_val = table[rows-1][j_ind]
			#self.printPointerMatrix(pointers)

			i = rows - 1
			j = j_ind
			score = table_min_val
			origin = (0,0)
			tuple = pointers[i][j]
			dist = (i, j) #goal
			align1 = ""
			align2 = ""

			#self.printTable(table)
			print(t_min_val_location)
			print("%s %s" %(i, j))
			print(dist)
			while(origin != dist):
				i_index = dist[0] - tuple[0]
				j_index = dist[1] - tuple[1]

		
				if(i_index == 1 and j_index == -1): # / banded = | unbanded
					ind_char_seq1 = dist[0] - 1
					align1 = seq1[ind_char_seq1] + align1
					align2 = "-" + align2
					
				if(i_index == 1 and j_index == 0): # | banded = \ unbanded
					ind_char_seq1 = dist[0] - 1
					ind_char_seq2 = (dist[0] + dist[1]) - 4 # (i + j) - 4
					align1 = seq1[ind_char_seq1] + align1
					align2 = seq2[ind_char_seq2] + align2
				
				if(i_index == 0 and j_index == 1): # __ banded and unbanded
					align1 = "-" + align1
					ind_char_seq2 = (dist[0] + dist[1]) - 4 # (i + j) - 4
					align2 = seq2[ind_char_seq2] + align2

				dist = tuple
				tuple = pointers[dist[0]][dist[1]]

			align1 = align1[:100]
			align2 = align2[:100]
			#print(align1)
			#print(align2)
		return score, align1, align2








	def printInfo(self, seq1, seq2, align_length, banded = False):
		if(banded):
			print("We are going to use banded alignment!")
			print("Sequence 1: " + seq1)
			print("Sequence 2: " + seq2)
			print("Alignment length: %s" % (align_length))
		else:
			print("we are using unrestricted alignment using Needleman-Wunsch Algorithm")
			print("Sequence 1: " + seq1)
			print("Sequence 2: " + seq2)
			print("Alignment length: %s" % (align_length))

	def printBandedInf(self, seq1, seq2, rows, cols):
		print("This is the inf for banded algorithm!")
		print("Sequence 1: " + seq1 + " with length of %s" % len(seq1))
		print("Sequence 2: " + seq2 + " with length of %s" % len(seq2))
		print("Table will have %s rows and %s cols" % (rows, cols))
		print("rows: %s , cols: %s" % (rows, cols))

	def printMatrix(self, matrix):
		print('\n'.join([' '.join(['{:4}'.format(item) for item in row]) for row in matrix]))
		print('\n')

	def printTable(self, matrix):
		print('\n'.join([' '.join(['{:4}'.format(item) for item in row]) for row in matrix]))
		print('\n')


	def printT(self, table):
		for i in range(len(table)):
			str = ' '
			for j in  range(len(table[i])):
				str = str + " " + repr(table[i][j])
			print(str)
		print()


	def printPointerMatrix(self, matrix):
		for i in range(len(matrix)):
			str = ''
			for j in range(len(matrix[i])):
				str = str + " " + self.convertTuple(matrix[i][j])
			print(str)
		print()
	def convertTuple(self, tup):
        # initialize an empty string
		str = "("
		index = 0
		for item in tup:
			if(index == len(tup)-1):
				str = str + (",%s" % item)
			else:
				str = str + ("%s" % item)
			index +=1
		str = str + ")"
		return str



#	x    0    p    o    l    y    n    o    m    i    a    l 
#	0    0    5   10   15   20   25   30   35   40   45   50
#	p    5   -3    2    7   12   17   22   27   32   37   42
#	o   10    2   -6   -1    4    9   14   19   24   29   34
#	l   15    7   -1   -9   -4    1    6   11   16   21   26
#	y   20   12    4   -4  -12   -7   -2    3    8   13   18
#	n   25   17    9    1   -7  -15  -10   -5    0    5   10
#	o   30   22   14    6   -2  -10  -18  -13   -8   -3    2
#	m   35   27   19   11    3   -5  -13  -21  -16  -11   -6
#	i   40   32   24   16    8    0   -8  -16  -24  -19  -14
#	a   45   37   29   21   13    5   -3  -11  -19  -27  -22
#	l   50   42   34   26   18   10    2   -6  -14  -22  -30


#	 	 0    e    x    p    o    n    e    n    t    i    a    l	
#	0    0    5   10   15   20   25   30   35   40   45   50   55
#	p    5    1    6    7   12   17   22   27   32   37   42   47
#	o   10    6    2    7    4    9   14   19   24   29   34   39
#	l   15   11    7    3    8    5   10   15   20   25   30   31
#	y   20   16   12    8    4    9    6   11   16   21   26   31
#	n   25   21   17   13    9    1    6    3    8   13   18   23
#	o   30   26   22   18   10    6    2    7    4    9   14   19
#	m   35   31   27   23   15   11    7    3    8    5   10   15
#	i   40   36   32   28   20   16   12    8    4    5    6   11
#	a   45   41   37   33   25   21   17   13    9    5    2    7
#	l   50   46   42   38   30   26   22   18   14   10    6   -1


#  Pointers Matrix for polynmial and exponential
#  (0,0) (0,0) (0,1) (0,2) (0,3) (0,4) (0,5) (0,6) (0,7) (0,8) (0,9) (0,10)
#  (0,0) (0,0) (1,1) (0,2) (1,3) (1,4) (1,5) (1,6) (1,7) (1,8) (1,9) (1,10)
#  (1,0) (1,1) (1,1) (2,2) (1,3) (2,4) (2,5) (2,6) (2,7) (2,8) (2,9) (2,10)
#  (2,0) (2,1) (2,2) (2,2) (3,3) (2,4) (3,5) (3,6) (3,7) (3,8) (3,9) (2,10)
#  (3,0) (3,1) (3,2) (3,3) (3,3) (4,4) (3,5) (4,6) (4,7) (4,8) (4,9) (4,10)
#  (4,0) (4,1) (4,2) (4,3) (4,4) (4,4) (5,5) (4,6) (5,7) (5,8) (5,9) (5,10)
#  (5,0) (5,1) (5,2) (5,3) (5,3) (5,5) (5,5) (6,6) (5,7) (6,8) (6,9) (6,10)
#  (6,0) (6,1) (6,2) (6,3) (6,4) (6,5) (6,6) (6,6) (7,7) (6,8) (7,9) (7,10)
#  (7,0) (7,1) (7,2) (7,3) (7,4) (7,5) (7,6) (7,7) (7,7) (7,8) (7,9) (8,10)
#  (8,0) (8,1) (8,2) (8,3) (8,4) (8,5) (8,6) (8,7) (8,8) (8,8) (8,9) (9,10)
#  (9,0) (9,1) (9,2) (9,3) (9,4) (9,5) (9,6) (9,7) (9,8) (9,9) (9,9) (9,10)




#			 -	 -	 -	 0 	 P 	 O	 L
#	   0  0 inf inf inf  0   5   10  15        Way to access char in seq1 is seq1(i - 1)
#	   1  P inf inf  5   *					   Way to access char in seq2 is seq2((i + j) - 4)
#	   2  O inf  10  
#	   3  L  15
#         Y
#    	  N
#  		  O
#		  M
#		  I
#		  A
#         L


#			 -	 -	 - 	 P 	 O	 L   Y
#	   0  P inf inf inf  0   5   10  15        Way to access char in seq1 is seq1(i)
#	   1  O inf inf  5   *					   Way to access char in seq2 is seq2((i + j) - 3)
#	   2  L inf  10  
#	   3  Y  15
#         N
#  		  O
#		  M
#		  I
#		  A
#         L


# rows = len(seq1)
# 			cols = bandWith
# 			self.printBandedInf(seq1, seq2, rows, cols)
# 			table = [[float("inf") for i in range(cols)] for j in range(rows)] # inf middle first row
# 			#self.printTable(table)
# 			self.printT(table)
# 			#range(startindex, stopindex (not inclusive), step)
# 			for j in range(3, -1, -1): # O(4)
# 				table[3-j][j] = (3-j)*5 # /

# 			for i in range(1): # O(1)
# 				for j in range(4): # O(4)
# 					table[i][j + 3] = j * 5 # __
			

# 			start = 4
			
# 			diagonal = float("inf")
# 			left = float("inf")
# 			top = float("inf")
# 			options = [] # left, top, diagonal 
# 			for i in range(1, rows): # O(n-1) = O(n)
# 				if(start > 0):
# 					start -= 1
# 				for j in range(start, cols, 1): # O(7) = O(1)
# 					index_of_seq2 = (i + j) - 3  #problem because 
# 					if(seq1[i] == seq2[index_of_seq2]):
# 						print("char in seq1 " + seq1[i] + " " "char in seq2 " + seq2[index_of_seq2])
# 						diagonal = table[i-1][j] - 3
# 					else:
# 						diagonal = table[i-1][j] + 1
					
# 					if(j != cols - 1): # right limit
# 						top = table[i-1][j+1] + 5

# 					if(j != 0): #left limit
# 						left = table[i][j-1] + 5
					
# 					options.append(left)
# 					options.append(top)
# 					options.append(diagonal)
# 					min_v = min(options)#O(3) = O(1)