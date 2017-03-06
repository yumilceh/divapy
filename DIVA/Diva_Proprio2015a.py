'''
Created on Feb 5, 2016
This sensorimor system defines the DIVA agent used for the CCIA 2015's paper
@author: Juan Manuel Acevedo Valle
'''

#import sys
#import wave
import subprocess as sp
import math
import numpy as np
import pymatlab as ml
import matplotlib.pyplot as plt
from matplotlib import animation 
from scipy.integrate import odeint
from scipy import linspace
from scipy.io.wavfile import write
#from matplotlib.pyplot import autoscale
#from matplotlib.animation import Animation
#from scipy.interpolate.interpolate_wrapper import block

class Diva_Proprio2015a:
    
    def __init__(self):
        motor_names=['M1','M2','M3','M4','M5','M6','M7','M8','M9','M10','M11','M12','M13','M14','M15','M16','M17','M18','M19','M20','M21','M22','M23','M24','M25','M26']
        sensor_names=['S1','S2','S3','S4','S5','S6']
        somato_names=['P1']
        n_motor=26
        n_sensor=6
        n_somato=1
        outputScale=[100.0,500.0,1500.0,3000.0];
        min_motor_values=np.array([-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-0.25,-0.25,-0.25]*2)
        max_motor_values=np.array([3.0,3,3,3,3,3,3,3,3,3,1,1,1]*2)
        
        min_sensor_values=np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        max_sensor_values=np.array([2.0, 2.0, 1.0, 2.0, 2.0, 1.0])
        
        ts=0.01;
        self.ts=ts;
        self.time=linspace(0, .8, int(.8/ts)+1)
        self.n_motor=n_motor
        self.n_sensor=n_sensor
        self.n_somato=n_somato
        self.motor_names=motor_names
        self.sensor_names=sensor_names
        self.somato_names=somato_names
        
        self.min_motor_values=min_motor_values
        self.max_motor_values=max_motor_values
        self.min_sensor_values=min_sensor_values
        self.max_sensor_values=max_sensor_values
        
        self.motor_command=[0.0] * n_motor
        self.sensorOutput=[0.0] * n_sensor
        self.sensor_goal=[0.0] * n_sensor
        self.somatoOutput=[0.0] * n_somato
        self.competence_result=0.0;
        self.matlabSession=ml.session_factory()        
        self.matlabSession.run('cd /home/yumilceh/eclipse_ws/Early_Development/SensorimotorExploration/Systems/DIVA/') #Path to DIVA functions
        self.matlabSession.putvalue('outputScale', outputScale)
        
    def setMotorCommand(self,motor_command):
        self.motor_command=motor_command    
        
    def getMotorDynamics(self,sound=0):
        if sound:
            ts=0.005;
        else:
            ts=self.ts
            
        durationM1=0.25
        durationM2=0.55
        nSamples=int(0.8/ts+1)
        nSamples1=int(durationM1/ts)+1
        nSamples2=int(durationM2/ts)+1
        y_neutral=[0.0]*13
        y_neutral[11]=-0.25
        y_neutral[12]=-0.25
        y0=[0.0]*26
        y0[:13]=y_neutral
        m1=self.motor_command[:13]
        t1=linspace(0.0,durationM1,nSamples1)
        artStates1=odeint(motorDynamics,y0,t1,args=(self,m1))
        t2=linspace(0.0,durationM2,nSamples2)
        m2=self.motor_command[13:]
        artStates2=odeint(motorDynamics,artStates1[-1,:],t2,args=(self,m2))
        if sound:
            return np.concatenate((artStates1,artStates2))
        else:
            self.artStates= np.zeros((nSamples, 26))
            self.artStates[:nSamples1,:]=artStates1
            self.artStates[nSamples1-1:,:]=artStates2
    
    def vocalize(self):
        ts=self.ts;
        perceptionWindowDuration=0.15;
        perceptionTime=linspace(ts,perceptionWindowDuration, int(perceptionWindowDuration/ts))
        nPerceptionSamples=(len(perceptionTime))
        self.auditoryResult=[0.0]*6
        proprioceptiveAv=[0.0]*2
        self.matlabSession.putvalue('artStates',self.artStates)
        #self.matlabSession.run('save artStates.mat artStates')
        self.matlabSession.run('mscript_Aud_Proprio')
        auditoryStates=self.matlabSession.getvalue('auditoryStates')
        self.auditoryStates=auditoryStates;
        minaf=self.matlabSession.getvalue('minaf')
        self.somatoOutput=minaf
        '''print('audStates')
        print(auditoryStates)
        print('minaf')
        print(minaf)'''
        
        #First perception time window
        for index in range(nPerceptionSamples):
            #print(index)
            #print(nPerceptionSamples)
            if (self.artStates[index+26,11]>0) and (self.artStates[index+26,12]>0) and (minaf[index+26]>0):
                self.auditoryResult[0]=self.auditoryResult[0]+auditoryStates[index+26,1]
                self.auditoryResult[1]=self.auditoryResult[1]+auditoryStates[index+26,2]
                self.auditoryResult[2]=self.auditoryResult[2]+1.0
            proprioceptiveAv[0]=proprioceptiveAv[0]+(minaf[index+26]/nPerceptionSamples)
        self.auditoryResult[0]=self.auditoryResult[0]/nPerceptionSamples
        self.auditoryResult[1]=self.auditoryResult[1]/nPerceptionSamples
        self.auditoryResult[2]=self.auditoryResult[2]/nPerceptionSamples   
            
        #Second perception time window
        for index in range(nPerceptionSamples):
            #print(index)
            if (self.artStates[index+66,11]>0) and (self.artStates[index+66,12]>0) and (minaf[index+66]>0):
                self.auditoryResult[3]=self.auditoryResult[3]+auditoryStates[index+66,1]
                self.auditoryResult[4]=self.auditoryResult[4]+auditoryStates[index+66,2]
                self.auditoryResult[5]=self.auditoryResult[5]+1.0
            proprioceptiveAv[1]=proprioceptiveAv[1]+(minaf[index+66]/nPerceptionSamples)
        self.auditoryResult[3]=self.auditoryResult[3]/nPerceptionSamples
        self.auditoryResult[4]=self.auditoryResult[4]/nPerceptionSamples
        self.auditoryResult[5]=self.auditoryResult[5]/nPerceptionSamples
                
        self.somatoOutput=0.0
        if((proprioceptiveAv[0]<0.0) or (proprioceptiveAv[1]<0.0)):
            self.somatoOutput=1.0
        self.sensorOutput=self.auditoryResult;     
            
    def executeMotorCommand(self):
        self.getMotorDynamics()
        self.vocalize()
    
    def plotArticulatoryEvolution(self,arts):
        for index in range(len(arts)):
            plt.plot(self.time,self.artStates[:,arts[index]-1])
            plt.hold(True)
        plt.show()
        
    def plotAuditoryOutput(self,audOut):
        for index in range(len(audOut)):
            plt.plot(self.time,self.auditoryStates[:,audOut[index]-1])
            plt.hold(True)
        plt.show()
        
    def plotSomatoOutput(self):
        plt.plot(self.time,self.somatoOutput)
        plt.show;
           
    def getVocaltractShape(self,artStates,returnShape=0):
        self.matlabSession.putvalue('artStates',artStates)
        self.matlabSession.run('[~, ~, outline] = diva_synth(artStates\', \'audsom\')')
        if returnShape:
            return self.matlabSession.getvalue('outline')
        else:
            self.vocalTractshape=self.matlabSession.getvalue('outline')
            
            
    def plotVocalTractShape(self,time,plot=1):
        self.getVocaltractShape(self.artStates)
        ts=self.ts
        index=np.round(time/ts)
        figVocalTract=plt.figure();
        axVocalTract=figVocalTract.add_subplot(111,autoscale_on=False, xlim=(-75,225),ylim=(-200,100))
        axVocalTract.plot(np.real(self.outline[:,index]), np.imag(self.outline[:,index]))
        figVocalTract.show()
        '''try:
            input('Press any key to continue...')
        except:
            pass'''
        
    def getVocalizationVideo(self,show=0):
        Writer = animation.writers['ffmpeg']
        writer = Writer(fps=1/0.005, metadata=dict(artist='Juan Manuel Acevedo Valle'))
        figVocalTract=plt.figure()
        
        artStates=self.getSoundWave(save=1,returnArtStates=1)
        outline=self.getVocaltractShape(artStates, returnShape=1)
        nSamples=artStates.shape[0]
        print(nSamples)
        sequence=[]
        for index in range(nSamples):
            sequence.append((plt.plot(np.real(outline[:,index]), np.imag(outline[:,index]))))
        im_ani = animation.ArtistAnimation(figVocalTract, sequence, interval=0.005, repeat=False,blit=True)
        im_ani.save('vt.mp4', writer=writer)
        command = ["ffmpeg",
                   '-i', 'vt.wav',
                   '-i','vt.mp4',
                   '-vcodec', 'copy', 'vtaudio.mp4']
        sp.call(command)
        if(show):
            figVocalTract.show();
            
        
               
    def getSoundWave(self, play=0,save=0,returnArtStates=0): #based on explauto
        soundArtStates=self.getMotorDynamics(sound=1)
        #print('ts=0.005')
        #print(soundArtStates.shape)
        #print('ts=0.01')
        #print(self.artStates.shape)
        self.matlabSession.putvalue('artStates',soundArtStates[:,0:13])
        #self.matlabSession.run('save artStates.mat artStates')
        self.matlabSession.run('soundWave = diva_synth(artStates\', \'sound\')')
        self.soundWave=self.matlabSession.getvalue('soundWave');
        if(play):
            self.playSoundWave()
        if(save):
            scaled = np.int16(self.soundWave/np.max(np.abs(self.soundWave)) * 32767)
            write('vt.wav', 11025, scaled)
        if(returnArtStates):
            return soundArtStates  

    def plotSoundWave(self):
        plt.plot(self.soundWave)
        plt.show();
    
    def playSoundWave(self): #keep in mind that DIVA works with ts=0.005
        import pyaudio 
        self.pa = pyaudio.PyAudio() #If pa and stream are not elements of the self object then sound does not play
        self.stream = self.pa.open(format=pyaudio.paFloat32,
                         channels=1,
                         rate=11025,
                         output=True)
        self.stream.start_stream()
        self.stream.write(self.soundWave.astype(np.float32).tostring())
    
    def releaseAudioDevice(self): #any sound in the buffer will be removed
        try:
            self.pa.terminate()
        except:
            pass
        
    def stop(self):
        del self.matlabSession

def motorDynamics(y,t,self,m):
    dumpingFactor=1.01
    w0=2*math.pi/0.8
    
    dy1=y[13]
    dy2=y[14]
    dy3=y[15]
    dy4=y[16]
    dy5=y[17]
    dy6=y[18]
    dy7=y[19]
    dy8=y[20]
    dy9=y[21]
    dy10=y[22]
    dy11=y[23]
    dy12=y[24]
    dy13=y[25]
        
    dy14=-2*dumpingFactor*w0*y[13]-(pow(w0,2))*y[0]+(pow(w0,2))*m[0]
    dy15=-2*dumpingFactor*w0*y[14]-(pow(w0,2))*y[1]+(pow(w0,2))*m[1]
    dy16=-2*dumpingFactor*w0*y[15]-(pow(w0,2))*y[2]+(pow(w0,2))*m[2]
    dy17=-2*dumpingFactor*w0*y[16]-(pow(w0,2))*y[3]+(pow(w0,2))*m[3]
    dy18=-2*dumpingFactor*w0*y[17]-(pow(w0,2))*y[4]+(pow(w0,2))*m[4]
    dy19=-2*dumpingFactor*w0*y[18]-(pow(w0,2))*y[5]+(pow(w0,2))*m[5]
    dy20=-2*dumpingFactor*w0*y[19]-(pow(w0,2))*y[6]+(pow(w0,2))*m[6]
    dy21=-2*dumpingFactor*w0*y[20]-(pow(w0,2))*y[7]+(pow(w0,2))*m[7]
    dy22=-2*dumpingFactor*w0*y[21]-(pow(w0,2))*y[8]+(pow(w0,2))*m[8]
    dy23=-2*dumpingFactor*w0*y[22]-(pow(w0,2))*y[9]+(pow(w0,2))*m[9]
    dy24=-2*dumpingFactor*w0*y[23]-(pow(w0,2))*y[10]+(pow(w0,2))*m[10]
    dy25=-2*dumpingFactor*w0*y[24]-(pow(w0,2))*y[11]+(pow(w0,2))*m[11]
    dy26=-2*dumpingFactor*w0*y[25]-(pow(w0,2))*y[12]+(pow(w0,2))*m[12]
    
    return [dy1,dy2,dy3,dy4,dy5,dy6,dy7,dy8,dy9,dy10,dy11,dy12,dy13,dy14,dy15,dy16,dy17,dy18,dy19,dy20,dy21,dy22,dy23,dy24,dy25,dy26]
    
