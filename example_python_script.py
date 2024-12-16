#!/usr/bin/env python3

import math,csv,glob
#for f1 in glob.glob('D:\\Cabo_mIF_QuPath\\QuPath\\S20-18383-A6'):
#F:\Cabo_Biopsy\C19-16520
for f1 in glob.glob('F:\\Cabo_Biopsy\\S21-20269'):
	print(f1)
	tcf_cd8=csv.reader(open(f1+'/tcf_cd8_cells.csv'))
	cd8=csv.reader(open(f1+'/cd8_cells.csv'))
	mhc=csv.reader(open(f1+'/mhc_cells.csv'))
	
	#example how to write density for each file
	#mhc_density=f1+'/mhcreal_density.csv'
	#mhc_writer=csv.writer(open(mhc_density,'w'))
	
	tcf_cd8_density=f1+'/tcf_cd8_density.csv'
	tcf_cd8_writer=csv.writer(open(tcf_cd8_density,'w'))
	
	cd8_density=f1+'/cd8_density.csv'
	cd8_writer=csv.writer(open(cd8_density,'w'))
	
	mhc_density=f1+'/mhc_density.csv'
	mhc_writer=csv.writer(open(mhc_density,'w'))
	
	
	line=0
	mhc=list(mhc)
	mhc_x=[]
	mhc_y=[]
	cd8=list(cd8)
	cd8_x=[]
	cd8_y=[]
	tcf_cd8=list(tcf_cd8)
	tcf_cd8_x=[]
	tcf_cd8_y=[]
	
	#2 and 3 are x and y columns (start counting from 0)
	for locations in mhc:
		#print(locations)
		if line>0 and locations:
			try:
				mhc_x.append(math.trunc(float(locations[6])))
				mhc_y.append(math.trunc(float(locations[7])))
			except:
				print('break22')
				#break
				#line+=1
		line+=1
		
	line=0
	print('test1')
	for locations in cd8:
		#print(locations)
		if line>0 and locations:
			try:
				cd8_x.append(math.trunc(float(locations[6])))
				cd8_y.append(math.trunc(float(locations[7])))
			except:
				print('break38')
				#break
				#line+=1
		line+=1
		
	line=0
	for locations in tcf_cd8:
		#print(locations)
		if line>0 and locations:
			try:
				tcf_cd8_x.append(math.trunc(float(locations[6])))
				tcf_cd8_y.append(math.trunc(float(locations[7])))
			except:
				print('break39')
				#break
				#line+=1
		line+=1
		
		
	mhc_l=zip(mhc_x,mhc_y)
	
	#for x in mhc_l:
	#	print(x)
	if len(mhc_x)<5:
		print('break29',len(mhc_x))
		break
	#Set window size and step length in pixels
	x_window=200
	y_window=200
	x_step=100
	y_step=100
	
	#Calculate # of windows across all samples
	total_x=tcf_cd8_x+mhc_x+cd8_x
	total_y=tcf_cd8_y+mhc_y+cd8_y
	
	x_dim=math.trunc(((max(total_x)-min(total_x))/x_step))+1
	y_dim=math.trunc(((max(total_y)-min(total_y))/y_step))+1
	print(x_dim,y_dim)
	
	
	mhc_list1=[]
	for y_position in range(min(total_y),max(total_y),y_step):
		for x_position in range(min(total_x),max(total_x),x_step):
			#print(y_position,x_position,'woot')
			count=0
			for x,y in zip(mhc_x,mhc_y):
				if int(x)>int(x_position) and int(x)<(int(x_position)+int(x_window)) and int(y)>int(y_position) and int(y)<(int(y_position)+int(y_window)):
					count +=1
					#print(count)
			mhc_list1.append(count)
			#print(count)
			#print(float(len(list1))/(x_dim*y_dim))
			
			
	cd8_l = zip(cd8_x,cd8_y)
	
	if len(cd8_x)<5:
		print('break114',len(cd8_x))
		break
	
	
	
	cd8_list1 = []
	for y_position in range(min(total_y),max(total_y),y_step):
		for x_position in range(min(total_x),max(total_x),x_step):
			#print(y_position,x_position,'woot')
			count=0
			for x,y in zip(cd8_x,cd8_y):
				if int(x)>int(x_position) and int(x)<(int(x_position)+int(x_window)) and int(y)>int(y_position) and int(y)<(int(y_position)+int(y_window)):
					count +=1
					#print(count)
			cd8_list1.append(count)
			#print(count)
			#print(float(len(list1))/(x_dim*y_dim))
			
	tcf_cd8_l = zip(tcf_cd8_x,tcf_cd8_y)
	
	if len(tcf_cd8_x)<5:
		print('break135',len(tcf_cd8_x))
		break
	
	
	tcf_cd8_list1 = []
	for y_position in range(min(total_y),max(total_y),y_step):
		for x_position in range(min(total_x),max(total_x),x_step):
			#print(y_position,x_position,'woot')
			count=0
			for x,y in zip(tcf_cd8_x,tcf_cd8_y):
				if int(x)>int(x_position) and int(x)<(int(x_position)+int(x_window)) and int(y)>int(y_position) and int(y)<(int(y_position)+int(y_window)):
					count +=1
					#print(count)
			tcf_cd8_list1.append(count)
			#print(count)
			#print(float(len(list1))/(x_dim*y_dim))
			
			
	i=0
	mhc_list2=[]
	
	#turn list into an x_dim x y_dim matrix
	while i<len(mhc_list1):
		mhc_list2.append(mhc_list1[i:i+x_dim])
		i+=x_dim
		
		#for z in mhc_list2:
		#mhc_writer.writerow(z)
		
	for z in mhc_list2:	
		mhc_writer.writerow(z)
		print('done')
		
	i = 0
	cd8_list2 = []
	
	# turn list into an x_dim x y_dim matrix
	while i < len(cd8_list1):
		cd8_list2.append(cd8_list1[i:i + x_dim])
		i+= x_dim
		
	for z in cd8_list2:
		cd8_writer.writerow(z)
		print('done')
		print('test2')
		
		
	i = 0
	tcf_cd8_list2 = []
	
	# turn list into an x_dim x y_dim matrix
	while i < len(tcf_cd8_list1):
		tcf_cd8_list2.append(tcf_cd8_list1[i:i + x_dim])
		i+= x_dim
		
	for z in tcf_cd8_list2:
		tcf_cd8_writer.writerow(z)
		print('done')