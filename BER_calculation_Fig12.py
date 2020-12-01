import pickle
import pdb
import hashlib
import scipy.io as sio

def hamming_distance(chaine1, chaine2):
    return sum(c1 != c2 for c1, c2 in zip(chaine1, chaine2))

def hamming_distance2(chaine1, chaine2):
    return len(list(filter(lambda x : ord(x[0])^ord(x[1]), zip(chaine1, chaine2))))

def rotate(strg, n):
    return strg[n:] + strg[:n]

def check_string(sentence, obj):
	try:
		return sentence.index(obj)
	except ValueError:
		return -1

def GenerateString(word,length):
	word_length = len(word)
	factor = length/word_length
	rem = length%word_length
	return word*factor + word[:rem]

def ByteCheck(sentence):
	length = len(sentence)
	byte_error = 10000000
	for x in range(0, 4):
		my_str = GenerateString(rotate("STOP",x),length)
		if my_str == sentence:
			return True
	return False

def BitCheck(sentence):
	length = len(sentence)
	bit_error = 10000000000
	sentence_binary = ''.join(format(ord(i), 'b') for i in sentence)
	for x in range(0, 4):
		my_str = GenerateString(rotate("STOP",x),length)
		my_str_binary = ''.join(format(ord(i), 'b') for i in my_str)
		current_bit_error = hamming_distance(my_str_binary,sentence_binary)
		if current_bit_error < bit_error:
			bit_error = current_bit_error
	return bit_error

nums = pickle.load( open( "iteration_nums_IROSFig12.pkl", "rb" ) )
data = pickle.load( open( "received_data_IROSFig12.pkl", "rb" ) )

LenDict = {}
ErrDict = {}
byte_error_array =[-1]*len(nums)
bit_error_array =[-1]*len(nums)
bit_error_rate = [2]*100
NumBitsArray = [0]*100
ErrBitsArray = [0]*100
length_array = [0]*len(nums)
for key in nums:
	length_array[key-1] = len(data[key]) 
	if ByteCheck(data[key]):
		byte_error_array[key-1] = 0
		bit_error_array[key-1] = 0 
	else:
		bit_error_array[key-1] = BitCheck(data[key])


	i = nums[key]
	if i in LenDict:
		LenDict[i] = LenDict[i] + length_array[key-1]
		ErrDict[i] = ErrDict[i] + bit_error_array[key-1]
	else:
		LenDict[i] = length_array[key-1]
		ErrDict[i] = bit_error_array[key-1]

	bit_error_rate[i] = 1.00000000*ErrDict[i]/(LenDict[i]*8)
	NumBitsArray[i] = LenDict[i]*8
	ErrBitsArray[i] = ErrDict[i]

sio.savemat("BitTransmissionDataIROS_Fig12.mat", {'NumBitsArray':NumBitsArray, 'ErrBitsArray':ErrBitsArray})

