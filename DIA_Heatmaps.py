# # Displaying the .mgf files as heatmaps using matplotlib

# I will be reading the data into a list of lists, rounding m/z to 3 decimals.
#
# Format is `[m/z value,[[spectrum 1,intensity]]],[m/z value,[[spectrum 2,intensity]]]`
#
# Then sort.
#
# Then merge into `[m/z value,[[spectrum 1,intensity],[spectrum 2,intensity],[etc]]]`
#
# Read those values into a numpy array.
#
# Plot using imshow.

import sys

#I have placed the entire script into a function so that I can run it in terminal
def main(data,pep): #pep is the peptide mass that signifies the swath
    import numpy as np
    import matplotlib.pyplot as plt
    import collections
    newData={}
    count=1
    size=0
    read=' '
    maxi=0


    '''
    While the read variable is still reading data in the file, get the spectrum number.
    If the spectrum has the correct pepmass, read in the data until the spectrum ends.
    Round the data to 3 decimals and put them into a dictionary with the keys being the m/z values.
    Print progress bar.
    '''

    dataset=open(data)

    print('Reading Data')

    while(''!=read):
        read=dataset.readline()
        if('Spectrum'+str(count) in read): #count is the spectrum number
            count+=1
            if(count%220==0): #progress bar
                print(int(count/220),end='...')
            if(count%2200==0):
                print()
        if('PEPMASS='+str(pep) in read): #checking for pepmass=pep
            size+=1
            while('END' not in read):
                read=dataset.readline()
                if('.' in read):
                    [x,y]=read.split(' ') #split the string of two numbers into x and y
                    x=np.round(float(x))
                    y=np.log2(np.round(float(y),3)) #round for computing purposes
                    if(y>maxi):
                        maxi=y
                    temp=[count-1,y] #create a list from the spectrum number and the intensity value (y)
                    newData.setdefault(x,[]).append(temp) #put that list to the appropriate m/z key (x)

    dataset.close()
    print()


    '''
    Run through the dictionary of lists and place the intensity values into an array.
    '''

    print('Sorting Data')

    newData=collections.OrderedDict(sorted(newData.items())) #this sorts the data by key value
    index=-1
    intensity=np.zeros([len(newData),size]) #create a new array of size m/z by number of spectra
    charge=np.zeros(len(newData))
    for key in newData.keys(): #for each m/z
        index+=1
        charge[index]=key
        for i in range(size): #for each spectrum
            try:
                intensity[index][(newData[key][i][0])]=((newData[key][i][1])/maxi) #put the intensity into the array
            except IndexError:
                break
    mini=charge[0]
    maxi=charge[-1]
    

    '''
    Plot
    '''

    print('Plotting Data')

    fig,ax=plt.subplots(figsize=(60,15))
    im=ax.imshow(intensity,cmap='Greys')

    y=ax.get_yticks()
    y=y/y[-1]*(maxi-mini)
    y+=mini
    y=np.round(y,0)

    ax.set_yticklabels(y)

    ax.set_title('Intensities of ions by m/z and spectrum number')
    fig.tight_layout()
    plt.savefig('Result_Graph.png')

# Command line usage
if __name__=='__main__': #Deja vu, I've just been in this place before, HIGHER ON THE STREET...
    main(sys.argv[1],sys.argv[2])
