import sys
import time
from datetime import datetime
import logging as log
import json
import glob
import random
from termcolor import colored

#import ROOT

#sys.path.append("/home/cbm/cbmsoft/emu_test_module_v2/python/smx_tester/")
import msts_defs as smc

class SmxCal(object):
    """class to provide calibration (trim) functionality to the smx class
       Implements methods for generating and checking calibrations
       Note: "set_trim" can be found in smx_oper  
    """
    def __init__(self):
        super().__init__()
        #self.trim_path = trim_path


    def get_trim_hash(self):
        """ calculate hash value of trim settings
        args: na
        returns: hash value
        """
        # ToDo: calculate hash also from trim FILE, not only from ASIC readback
        trim_l = [[disc for disc in range(31)] for ch in range(128)]
        for ch in range(128):
            for disc in range(31):
                trim_l[ch][disc] = self.read(ch, 2 * disc + 1) & 0xff
        trim_hash =  hash(tuple(map(tuple, trim_l)))     # ToDo: investigate other hashes from hashlib
        return trim_hash

    def get_trim_select_hash(self) -> int:
        """ calculate hash value of a SUBSET of trim settings; includes FAST disc trims
        args: na
        returns: hash value
        """
        disc_set = [0, 30, 33]
        trim_l = [[disc for disc in range(34)] for ch in range(128)]
        for ch in range(128):
            for disc in disc_set:
                trim_l[ch][disc] = self.read(ch, 2 * disc + 1) & 0xff
        trim_hash =  hash(tuple(map(tuple, trim_l)))     # ToDo: investigate other hashes from hashlib
        return trim_hash
        
    def scan_scurve(self, vp_min: int = 10, vp_max: int = 260, vp_step: int = 10, n_pulse: int = 30,
                    filename_scan: str = "scurves/scurves_"):
        """ perform s-curve scan
            copies and quickly adjusted to cri python sw
        args:    see above (TODO)
        returns: na
        saves raw scurve counter values to file
        """
        sn = 99
        amp_cal_min = 30
        amp_cal_max = 240
        amp_cal_fast = 40
        now = datetime.now()
        date =  hex(self.a_id)+now.strftime("-%y%m%d_%H%M%S")
        filename_scan = ( filename_scan + date
                          + "_sn_" + "{0:0>{1}}".format( sn, 3)
                          + "_asic_addr_" + str( self.hw_addr )
                          + "_" + "{0:0>{1}}".format(self.disc_th2_glob, 3)
                          + "_" + "{0:0>{1}}".format(self.vref_n, 3)
                          + "_" + "{0:0>{1}}".format(self.vref_p, 3)
                          + "_" + "{0:0>{1}}".format(self.vref_t, 3)
                          + "_vp_" + "{0:0>{1}}".format(amp_cal_min, 3)
                          + "_" + "{0:0>{1}}".format(int(amp_cal_max), 3)
                          + "_" + "{0:0>{1}}".format(amp_cal_fast, 3)    )

        # set polarity;  looks a bit self-consistent, but makes
        #   sure that configuration polarity is really written to ASIC
        self.conf_func("chrg_type", self.chrg_type)
        vp_nb = int( (vp_max - vp_min) / vp_step )
        cnt = [[[0 for vp in range(vp_nb)] for d in range(31)] for ch in range(128)]
        for vp_idx, vp in enumerate(range(vp_min, vp_max, vp_step)):
            # loop over groups
            for grp in range(0, 4):
                self.write(130, 5, grp)
                self.reset_adc_counter()
                self.gen_test_pulses(n_pulse, vp, grp)
                # readback counter values
                for d in range(0, 31):
                    n = 2 * d
                    for ch in range(grp, 128, 4):
                        cnt[ch][d][vp_idx] = self.read(ch, n)

        # write scan results in a format compatible with previous SW
        # ---------------  Polarity selection --------------
        if (self.chrg_type == 1):
            filename_scan = filename_scan + "_holes.txt"
        elif (self.chrg_type == 0):
            filename_scan = filename_scan + "_elect.txt"
        else:
            log.error("Undefined charge type!! Exiting")
            exit()
        myfile = open(filename_scan, "w+")
        # Loop on Channels
        for ch in range( 0, 128):
            myfile.write( "\n" )
            myfile.write( " ------------------------------------------ S-curves in channel: " + '{:3d}'.format(ch)
                          + "- --------------------------------------\n"  )
            myfile.write( "\n" )

            log_line = "                "
            for d in range(0, 31):
                log_line = log_line + '{:5d}'.format( d )
            myfile.write( log_line )
            myfile.write( "\n" )

            vp_idx = 0
            for vp in range(vp_min, vp_max, vp_step):
                log_line = "Volt_pulse " + '{:3d}'.format(vp) + ": "
                for d in range(0, 31):
                    log_line = log_line + '{:5d}'.format( cnt[ch][d][vp_idx])
                myfile.write( log_line )
                myfile.write( "\n" )
                vp_idx += 1
            myfile.write( " -----------------------------------------------------------------------"
                          + " --------------------------------------" )
            myfile.write( "\n" )
        myfile.close()
   
    def get_trim_adc_SA(self, pol, trim_final: list, loop_number: int = 40, amp_cal_min: int = 30, amp_cal_max: int = 247, much_mode_on: int = 0):

        loop_max = loop_number
        loop_max_fine = loop_number
        tr_coarse_step = 5
        tr_coarse_range = 1
        itr_coarse = 12
        rng = 25
        itr_fine = 20
        tr_i_fine_max = itr_fine
        n_length = itr_coarse
        width = rng
        tr_fine_offset = -15
        tr_max= 63
        
        # ------------------ get_trim_settings ------------------------------------------------------
        grp_min = 0
        grp_max = 4
        #channel range
        ch_min = 0
        ch_max = 128
        #discriminator
        d_min = 0
        d_max = 31

        #tr_i_fine_max =tr_coarse_step*tr_coarse_range-tr_fine_offset
        tr_fine_step = 1

        #----------------pol------------------------------
        if (pol == 1):
            self.write(130,2,163)   #  163 holes     131 electrons
            log.info(" ")
            log.info("____HOLES MODE____")
            log.info(" ")
        if (pol == 0):
            self.write(130,2,131)   #  163 holes     131 electrons
            log.info(" ")
            log.info("____ELECTRONS MODE____")
            log.info(" ")

        #----------MUCH-----------------
        #Read 130,1
        old_pulser_mode = ( self.read(130, 1) & 0xBF )
        if (much_mode_on == 1):
            log.info("____MUCH MODE____")
            log.info(" ")
            #Force MUCH mode for all channels
            for ch in range( 0, 130 ):
                self.write(ch, 65, 228 )
            #Set 130,1 for pulser in MUCH mode
            self.write( 130, 1, old_pulser_mode | 0x40 )
            
        # -----------------------------------------------------------------
        #print " "
        #print ".......................... Getting ADC Trim Values ................................"
        #print " "

        #Counters array for trim
        vpset = [0 for d in range (d_min,d_max)]   # array of pulse heights for each discriminator
        cnt = [[0 for d in range(d_max)] for ch in range(ch_max)]
        aux = [[100 for d in range (d_max)] for ch in range (ch_max)]
        trim_before = [[128 for d in range(d_max)] for ch in range(ch_max)]
        #  trim_low = [[[0 for tr in range(tr_i_fine_max+1)] for d in range(d_max)] for ch in range(ch_max)]        # lower trim values before switching
        #  cnts_aux = [[[0 for tr in range(tr_i_fine_max+1)] for d in range(d_max)] for ch in range(ch_max)]
        value1 = [[0 for d in range (d_max)] for ch in range (ch_max)]
        fcnt = [[[0 for tr in range(tr_max)] for d in range(d_max)] for ch in range(ch_max)]
        hh_cnt = [[0 for d in range(d_max)] for ch in range(ch_max)]
        #time_ADC_COARSE = 0
        #time_ADC_FINE   = 0
        vp_d=((amp_cal_max-amp_cal_min)/(d_max-d_min))# + 0.5 )    # + 0.5 Why ?     # Step in the pulse amplitud/trim
        for d in range (d_min, d_max):
            vpset[d]= int(amp_cal_min + (d-d_min)*vp_d)
        # Read slow shaper configuration
        shslowfs = ( self.read(130, 5) & 0xC ) >> 2
        # Loop over discriminators
        for d in range (d_min, d_max):
            disc = 61 - 2*d
            count = 61 - 2*d -1
            #vp = int(vpset[d])
            if (vpset[d]<amp_cal_min or vpset[d]>amp_cal_max):
                print ("NOTE: Pulse amplitude should be in the range ....")
            self.write(130, 4, vpset[d])
            print ("Calibration Pulse Amplitude set to: ", vpset[d])
            # TIME ADC COARSE
            #time_start_ADC_COARSE = time.time()
            for n in range (1, n_length):
                self.write(192, 2, 32)
                self.write(192, 2, 0)
                for grp in range (grp_min, grp_max):
                    grp_shslow = ((shslowfs & 0x3) << 2 | (grp & 0x3))
                    self.write(130, 5, grp_shslow) 
                    for loop in range (0, loop_max):
                        self.write(130, 11, 128)
                        self.write(130, 11, 0)
                for ch in range (ch_min, ch_max):
                    cnt[ch][d] = self.read(ch, count) & 0xfff
                    if (cnt[ch][d] > loop_max/2.):
                        trim_before[ch][d] = trim_before[ch][d] - int(width/n)
                        #trim_coarse_low[ch][d] = trim_before[ch][d]
                        trim_final[ch][d] = trim_before[ch][d]
                        #aux[ch][d] = abs(cnt[ch][d] - loop_max/2.)
                        self.write(ch, disc, trim_before[ch][d])
                    elif (cnt[ch][d] < loop_max/2.):
                        trim_before[ch][d] = trim_before[ch][d] + int(width/n)
                        self.write(ch, disc, trim_before[ch][d])
                        #trim_coarse_low[ch][d] = trim_before[ch][d]
                        trim_final[ch][d] = trim_before[ch][d]
                        #aux[ch][d] = abs(cnt[ch][d] - loop_max/2.)
                    else:
                        self.write(ch, disc, trim_before[ch][d])
                        #trim_coarse_low[ch][d] = trim_before[ch][d]
                        trim_final[ch][d] = trim_before[ch][d]
                        #aux[ch][d] = abs(cnt[ch][d] - loop_max/2.)
    
            #fine loop over tr#
            for itr in range(0,tr_i_fine_max+1):
                log.info("fine_disc %3d itr %3d/%3d", d, itr, tr_i_fine_max)
                ## Loop over groups
                for grp in range (grp_min, grp_max):
                    grp_shslow = ((shslowfs & 0x3)<<2 | (grp & 0x3)) # Included the slow shaper configuration (90,130,220,280ns)
                    self.write(130, 5, grp_shslow)
                    
                    # Loop over channels in group to set trim values
                    for ch in range(grp,ch_max,4):
                        tr = trim_before[ch][d] + itr*tr_fine_step - 5
                        self.write(ch,disc,tr)


                    ## Reset ADC counters
                    self.write(192,2,32)
                    #time.sleep(0.00001)
                    self.write(192,2,0)

                    ## multiple count loops
                    for loop in range(0,loop_max_fine):
                        self.write(130,11,128)    # sending trigger pulses ()
                        #time.sleep(0.004)
                        self.write(130,11,0)
                                
                    ## loop to read counter
                    #cnt_val = 0
                    for ch in range(grp,ch_max, 4):
                        fcnt[ch][d][itr] = self.read(ch,count) 
                    #log.info("\n"

            #Smooth fcnt -> This is not done here
            asum = 0
            isum = 0
            avg_max = 0.000
            avg_max_range = 5
            #logfile.write("the Half height counts condition")
            for ch in range(ch_min,ch_max):
                avg_max = 0
                for itr in range (tr_i_fine_max-avg_max_range,tr_i_fine_max):
                    #avg_cnt[ch][d][itr] = fcnt[ch][d][itr]
                    #if (itr>=(tr_i_fine_max-avg_max_range) and (itr<(tr_i_fine_max))):
                    if (abs(fcnt[ch][d][itr]-loop_max_fine)>0.1*loop_max_fine):
                        fcnt[ch][d][itr] = loop_max_fine
                    avg_max += fcnt[ch][d][itr]      # Max value condition
                hh_cnt[ch][d] = (avg_max/avg_max_range/2.)   # Condition for half height count

                # determining switching point for trim_final[ch][d]
                #logfile.write("Determining switching point for final trim")
                #for ch in range(ch_min,ch_max):
                trim_final[ch][d] = -1
                find_flag = 0
                y_min = 0
                y_max = 0
                for itr in range(0,tr_i_fine_max):
                    if (itr > 0 and itr<(tr_i_fine_max-1) and find_flag ==0):
                        if ((fcnt[ch][d][itr]<=hh_cnt[ch][d]) and (fcnt[ch][d][itr+1]>=hh_cnt[ch][d])):
                            y_max = fcnt[ch][d][itr+1] - hh_cnt[ch][d]
                            y_min = hh_cnt[ch][d] - fcnt[ch][d][itr]
                            trim_final[ch][d] = trim_before[ch][d] - 5 + ((itr*tr_fine_step) if (y_min<=y_max) else (itr+1)*tr_fine_step)
                            find_flag =1

        # restore original setting
        if (much_mode_on == 1):
            self.write( 130, 1, old_pulser_mode )
                    
        return 0

    def get_trim_adc(self, pol, trim_final: list, npulses: int = 40, amp_cal_min: int = 30, amp_cal_max: int = 247, vref_t: int = 118, much_mode_on: int = 0):

        # ------------------ get_trim_settings ------------------------------------------------------                        
        grp_min = 0
        grp_max = 4
        # channel range                                                                                                      
        ch_min = 0
        ch_max = 128
        # discrim. range. Here fast disc is not included                                                                     
        d_min = 0
        d_max = 31

        amp_cal_min = amp_cal_min
        amp_cal_max = amp_cal_max
        tr_min = 30
        tr_max = 220
        tr_coarse_step = 5
        tr_coarse_range = 1
        tr_coarse_offset = -5
        tr_fine_offset = -20

        loop_max = npulses
        loop_max_fine = loop_max
        vref_t = vref_t
        # ................. ADC thr settings (0...137,d_min...d_max,thr) ..............                                      
        itr = 0
        tr_i_fine_max =tr_coarse_step*tr_coarse_range-tr_fine_offset+15
        tr_fine_step = 1
        
        # calibration range                                                                                                  
        cnt_min_coarse = int(loop_max*0.30)
        cnt_max_coarse = int(loop_max*0.40)
        
        # ------------------------ pol -------------------------------------                                                 
        if (pol == 1):
          self.write_check(130,2,163)   #  163 holes     131 electrons                                          
          log.info(" ")
          log.info("____HOLES MODE____")
          log.info(" ")
        if (pol == 0):
          self.write(130,2,131)   #  163 holes     131 electrons                                                
          log.info(" ")
          log.info("____ELECTRONS MODE____")
          log.info(" ")

          self.write(130,18, vref_t)
        # -----------------------------------------------------------------                                                  
        # ------------------------ MUCH -----------------------------------                                                  

        # Read 130,1 configuration                                                                                           
        old_pulser_mode = ( self.read(130, 1) & 0xBF )
        if (much_mode_on == 1):
            log.info("____MUCH MODE____")
            log.info(" ")
            # Force MUCH mode for all channels                                                                                 
            for ch in range( 0, 130 ):
                self.write_check(ch, 65, 228 )
                # Set 130,1 for pulser in MUCH mode                                                                                
            self.write_check( 130, 1, old_pulser_mode | 0x40 )
                
        # -----------------------------------------------------------------

        log.info(" ")
        log.info(".......................... Getting ADC Trim Values ................................")
        log.info(" ")
        # counters array for trim                                                                                                                                                                                        
        vpset = [0 for d in range (d_min,d_max)]                                                                                                  # array of pulse heights for each discriminator                        
        vcnt = [[[0 for tr in range(tr_max)] for d in range(d_max)] for ch in range(ch_max)]    # array for discriminator counters coarse                                                                                
        fcnt = [[[0 for tr in range(tr_max)] for d in range(d_max)] for ch in range(ch_max)]          # array for discriminator counters fine                                                                            
        avg_cnt = [[[0 for tr in range(tr_max)] for d in range(d_max)] for ch in range(ch_max)] # smothed counter values from fine scan                                                                                  
        hh_cnt = [[0 for d in range(d_max)] for ch in range(ch_max)]                                                          # half maximun count values                                                                
        trim_coarse_low = [[0 for d in range(d_max)] for ch in range(ch_max)]                                           # lower trim values before switching                                                             
        
        # setting vpset = [32]                                                                                                                                                                                           
        # Implemented linearily -> could it be not linearly to implemment nonlinearADC characteristics                                                                                                                   
        vp_d=((amp_cal_max-amp_cal_min)/(d_max-d_min))# + 0.5 )    # + 0.5 Why ?     # Step in the pulse amplitud/trim                                                                                                   
        log.info("vp_d : %3f", float(vp_d) )
        
        for d in range (d_min, d_max):
            vpset[d]= int(amp_cal_min + (d-d_min)*vp_d)
            log.info("vpset: %3f", float(vpset[d]))

        # Read slow shaper configuration                                                                                                                                                                                 
        shslowfs = ( self.read(130, 5) & 0xC ) >> 2
        d_counter = 0
        # Loop over discriminators
        for d in range (d_min, d_max):
            disc = 61 - 2*d
            count = 61 - 2*d -1
            vp = int(vpset[d])
            if (vp<amp_cal_min or vp>amp_cal_max):
                log.info("NOTE: Pulse amplitude should be in the range ....")
                self.write(130, 4, vp)
                log.info("Calibration Pulse Amplitude set to: %3d", vp)
                itr = 0
                #Coarse loop over trim                                                                                                                                                                                         
            for tr in range (tr_min,tr_max,tr_coarse_step):
                #log.info("Discriminator number:   " , d, "  Set_trim:   ", tr, "\n")                                                                                                                                        
                #if tr<0 or tr>255:                                                                                                                                                                                          
                #log.info("NOTE: Trim out of range, should be between (0-255) ....", tr, "\n")                                                                                                                             
                # Loop over groups                                                                                                                                                                                         
                for grp in range (grp_min, grp_max):
                    grp_shslow = ((shslowfs & 0x3)<<2 | (grp & 0x3)) # Included the slow shaper configuration (90,130,220,280ns)                                                                                               
                    self.write(130, 5, grp_shslow)
                    #log.info(" Selected Shaping time and group: " ,format(bin(sts.read(130,5))), "\n")                                                                                                                        
                    # Loop over the channels of the group                                                                                                                                                                      
                    for ch in range(grp, ch_max, 4):
                        self.write(ch, disc, tr)
                    # Reset ADC counters                                                                                                                                                                                     
                    self.write(192,2,32)
                    #time.sleep(0.00001)                                                                                                                                                                                     
                    self.write(192,2,0)
                    # multiple count loops                                                                                                                                                                                   
                    for loop in range (0,loop_max):
                        #Sending the signal to trigger pulses.                                                                                                                                                                   
                        self.write(130,11,128)
                        #time.sleep(0.004)                                                                                                                                                                                       
                        self.write(130,11,0)
                    # loop to read counter                                                                                                                                                                                     
                    cnt_val = 0
                    for ch in range(grp,ch_max, 4):
                        vcnt[ch][d][itr]= self.read(ch,count)
                        #vcnt[ch][d][itr]= (self.read(ch,count) & 0x3FF)      # this is for the stacked bit ASIC                                                                                                    
                itr +=1  # here ends the loop over trim

            # Finding coarse switching point from vcnt[ch][d][tr]                                                                                                                                                          
            #log.info("............Coarse switching..........")                                                                                                                                                            
            #logfile.write(" ............Coarse switching..........\n")                                                                                                                                                    
            for ch in range(ch_min,ch_max):
                trim_coarse_low[ch][d] = 0                                  
                #trim_coarse_low[ch][d] = -1
                itr = 0
                #print "ch", '{:4d}'.format(ch),                                                                                           
                coarse_flag =0
                y_coarse_min =0
                y_coarse_max =0
                for tr in range (tr_min,tr_max,tr_coarse_step):   # the loop should include the tr_max
                    if (itr >= tr_coarse_range and coarse_flag == 0):
                        if ((vcnt[ch][d][itr]>=cnt_max_coarse) and vcnt[ch][d][itr-tr_coarse_range]<=cnt_max_coarse):
                            y_coarse_min = cnt_max_coarse-vcnt[ch][d][itr-tr_coarse_range]
                            y_coarse_max = vcnt[ch][d][itr] - cnt_max_coarse
                            trim_coarse_low[ch][d] = (tr if (y_coarse_min<=y_coarse_max) else (tr+tr_coarse_step)) - (tr_coarse_range*tr_coarse_step - tr_coarse_offset) # Searching range tr-25(35)
                            if (trim_coarse_low[ch][d]<0):
                                trim_coarse_low[ch][d]=0
                            coarse_flag =1
                    itr +=1
            log.info("")

            #fine loop over tr                                                                                      
            for  itr in range(0,tr_i_fine_max+1):
                log.info("fine_disc %3d itr %3d/%3d", d, itr, tr_i_fine_max)
                for grp in range (grp_min, grp_max):
                    grp_shslow = ((shslowfs & 0x3)<<2 | (grp & 0x3)) # Included the slow shaper configuration (90,130,220,280ns)
                    self.write(130, 5, grp_shslow)
                    #log.info(" Selected Shaping time and group: " ,format(bin(sts.read(130,5))), "\n"                  
                    # Loop over channels in group to set trim values                                                    
                    for ch in range(grp,ch_max,4):
                        tr = trim_coarse_low[ch][d] + itr*tr_fine_step
                        if (tr == 0):
                            tr = tr_min
                        log.info("Ch: {}  Disc: {} Thr: {}".format(ch, disc, tr))
                        self.write(ch,disc,tr)
                    ## Reset ADC counters                                                                               
                    self.write(192,2,32)
                    #time.sleep(0.00001)                                                                                
                    self.write(192,2,0)
                    ## multiple count loops                                                                             
                    for loop in range(0,loop_max_fine):
                        self.write(130,11,128)    # sending trigger pulses ()                                
                        #time.sleep(0.004)                                                                                
                        self.write(130,11,0)
                    ## loop to read counter                                                                             
                    #cnt_val = 0                                                                                        
                    for ch in range(grp,ch_max, 4):
                        fcnt[ch][d][itr]= self.read(ch,count)
                        #log.info("\n"                                     

            # Smooth fcnt -> This is not done here                                                                  
            asum = 0
            isum = 0
            avg_max = 0.000
            avg_max_range = 5
            #logfile.write("the Half height counts condition")                                                      
            for ch in range(ch_min,ch_max):
                avg_max = 0
                for itr in range (tr_i_fine_max-avg_max_range,tr_i_fine_max):
                    #avg_cnt[ch][d][itr] = fcnt[ch][d][itr]                                                             
                    #if (itr>=(tr_i_fine_max-avg_max_range) and (itr<(tr_i_fine_max))):                                 
                    if (abs(fcnt[ch][d][itr]-loop_max_fine)>0.1*loop_max_fine):
                        fcnt[ch][d][itr] = loop_max_fine
                    avg_max += fcnt[ch][d][itr]      # Max value condition                                            
                hh_cnt[ch][d] = (avg_max/avg_max_range/2.)   # Condition for half height count                        
                # determining switching point for trim_final[ch][d]                                                   
                #logfile.write("Determining switching point for final trim")                                          
                #for ch in range(ch_min,ch_max):
                trim_final[ch][d] = -1
                find_flag = 0
                y_min = 0
                y_max = 0
                for itr in range(0,tr_i_fine_max):
                    #print '{:4d}'.format(avg_cnt[ch][d][itr]),                                                         
                    if (itr > 0 and itr<(tr_i_fine_max-1) and find_flag ==0):
                        if ((fcnt[ch][d][itr]<=hh_cnt[ch][d]) and (fcnt[ch][d][itr+1]>=hh_cnt[ch][d])):
                            #print "   |   "                                                                                
                            y_max = fcnt[ch][d][itr+1] - hh_cnt[ch][d]
                            y_min = hh_cnt[ch][d] - fcnt[ch][d][itr]
                            trim_final[ch][d] = trim_coarse_low[ch][d] + ((itr*tr_fine_step) if (y_min<=y_max) else (itr+1)*tr_fine_step)
                            find_flag =1

        # restore original setting                                                                                                              
        if (much_mode_on == 1):
            self.write_check( 130, 1, old_pulser_mode )
            
        return 0
                        
                
    def get_trim_fast(self, pol, trim_final: list, npulses: int = 100, amp_cal_fast: int = 50, much_mode_on: int = 0 ):
        """ Determine full set of fast trim values
        legacy implementation from old SW. to be reviewed
        args: pls check code
        """
        grp_min = 0
        grp_max = 4

        ch_min = 0
        ch_max = 128

        thr_min= 0
        thr_max = 63
        thr_step = 1
        npulses = npulses
        # ------------------------ pol -------------------------------------

        if (pol == 1):
            # self.write(130,2,163)   #  163 holes     131 electrons
            log.info(" ")
            log.info("____HOLES MODE____")
            log.info(" ")
        if (pol == 0):
            # self.write(130,2,131)   #  163 holes     131 electrons
            log.info(" ")
            log.info("____ELECTRONS MODE____")
            log.info(" ")
        # -----------------------------------------------------------------
        # ------------------------ MUCH -----------------------------------

        # Read 130,1 configuration
        old_pulser_mode = ( self.read(130, 1) & 0xBF )
        if (much_mode_on == 1):
            log.info("____MUCH MODE____")
            log.info(" ")
            # Force MUCH mode for all channels
            for ch in range( 0, 130 ):
                self.write(ch, 65, 228 )
                # Set 130,1 for pulser in MUCH mode
            self.write( 130, 1, old_pulser_mode | 0x40 )

        cnt = [[0 for thr in range(thr_max)] for ch in range(ch_max)]
        avg = [[0 for thr in range(thr_max)] for ch in range(ch_max)]
        hh  = [0 for ch in range(ch_max)]

        thr_i = 0

        self.write(130, 4, amp_cal_fast)

        for thr in range(thr_min, thr_max):
            log.info("disc_threshold: %3d", thr )
            for grp in range(grp_min, grp_max):
                self.write(130, 5, grp)
                log.info("group: %3d", grp )

                # setting disc. threshold
                for ch in range(grp, ch_max, 4):
                    self.write(ch, 67, thr)

                # reseting counters
                self.write(192, 2, 32)
                self.write(192, 2, 0)

                # generating npulses
                for n in range(npulses):
                    self.write(130, 11, 128)
                    self.write(130, 11, 0)

                # reading ch counters
                for ch in range(grp, ch_max, 4):
                    cnt[ch][thr_i] = self.read(ch, 62)
                    log.info("FAST Calib:   ch: %3d %3d", ch, cnt[ch][thr_i])
                #      sys.stdout.flush()
                log.info(" ")
            thr_i += 1
            log.info(" ")

        thr_i = 0
        # thr_val_t = int(thr_max/thr_step)
        # ---------------------- Finding trim values -----------------------
        # -------------- Smoothing the  scurve and finding S-curves HH -----------------
        for ch in range(ch_min, ch_max):
            thr_i = 0
            isum = 0
            asum = 0
            avg_max = 0
            avg_max_range = 5
            # print "ch: ", ch, " ",
            for thr in range(thr_min, thr_max, thr_step):
                if (thr <= thr_step):
                    avg[ch][thr_i] = 0
                elif (thr >= thr_max-thr_step):
                    avg[ch][thr_i] = 0
                else:
                    isum = cnt[ch][thr_i-1]+cnt[ch][thr_i]+cnt[ch][thr_i+1]
                    asum = float(isum)
                    avg[ch][thr_i] = int(asum/3)

                if (thr >= thr_max - 6*thr_step and thr < thr_max-thr_step):
                    if (abs(cnt[ch][thr_i]-npulses) > 0.1*npulses):
                        cnt[ch][thr_i] = npulses
                    avg_max += cnt[ch][thr_i]
                thr_i += 1

            hh[ch] = int(int(avg_max/avg_max_range)/2)
            # print '{:3d}'.format(hh[ch])

        # print " "
        # print "------------ Final trim values ------------ "
        # print " "
        # print type(trim_final)

        for ch in range(ch_min, ch_max):
            thr_i = 0
            trim_final[ch][31] = -1
            find_flag = 0
            y_min = 0
            y_max = 0
            log_str = "CH " + str(ch).rjust(3)

            for thr in range(thr_min, thr_max, thr_step):
                if thr_i % 2 == 0:
                    log_str = log_str + "  " + str(cnt[ch][thr_i]).rjust(4)
                if (thr > 0 and thr < thr_max-thr_step and find_flag == 0):
                    # if (avg[ch][thr_i]<=hh[ch] and avg[ch][thr_i-1]<hh[ch] and avg[ch][thr_i+1]>= hh[ch]):
                    if (cnt[ch][thr_i] <= hh[ch] and cnt[ch][thr_i+1] >= hh[ch]):
                        #          print  " | "
                        y_max = cnt[ch][thr_i+1]-hh[ch]
                        y_min = hh[ch]-cnt[ch][thr_i]
                        trim_final[ch][31] = (thr_i*thr_step) if (y_min <= y_max) else ((thr_i+1)*thr_step)
                        find_flag = 1
                thr_i += 1
            # print "\n"
            log.info("%s", log_str)

        # restore original setting
        if (much_mode_on == 1):
            self.write( 130, 1, old_pulser_mode )

        return 0

    def write_trim_file(self, filename_trim: str, pol: int, trim_final: list, amp_cal_min: int = 40, amp_cal_max: int = 226, amp_cal_fast: int = 50, much_mode_on: int = 0):
        """ Write full set of trim values to file
        legacy implementation from old SW. to be reviewed
        args: polarity ( 0: electrons   1: holes )
              filename
              set of trims
              much_mode (default 0 - disabled)
        """

        pol = (self.read(130, 2) >> 5) & 0x1
        if pol == 0:
            polstr = "_elect"
        elif pol == 1:
            polstr = "_holes"
        else:
            polstr = "_na"

        # Much mode: calibration made for Much
        # ------------------------ MUCH -----------------------------------
        if (much_mode_on == 1):
            log.info(" ")
            log.info("____WRITING TRIMING FILE for MUCH ____")
            log.info(" ")
            filename_trim = filename_trim + "_much"
        # -----------------------------------------------------------------
        # Reading ASIC unique E-fused ID
        id_dev = self.read_efuse_str()
        trimfile = open(filename_trim + '.txt', "w+")
        # 11 lines of settings and calibration information before the trim
        trimfile.write("# ID " + str(id_dev) + "\n")
        trimfile.write("# Tag " + "_trim_" + "\n")
        trimfile.write("# Date " + datetime.now().strftime("%y%m%d_%H:%M:%S")+"\n")
        trimfile.write("# Vref_p " + str(self.read(130, 9) & 0xff) + "\n")
        trimfile.write("# Vref_n " + str(self.read(130, 8) & 0xff) + "\n")
        trimfile.write("# Thr2_glb " + str(self.read(130,7) & 0xff) + "\n")
        trimfile.write("# Vref_t "+str(self.read(130, 18) & 0xff) + "\n")
        trimfile.write("# Vref_t_range " + str((self.read(130, 18) & 0xc0) >> 6) + "\n")
        trimfile.write("# Pol " + str(pol) + "\n")
        trimfile.write("# ADC range " + str(amp_cal_min) + "-" + str(amp_cal_max) + "\n")
        trimfile.write("# FAST disc " + str(amp_cal_fast) + "\n")
        
        log.info("")
        log.info("")
        # writing trim values on file
        for ch in range(128):
            trimfile.write("ch:")
            trimfile.write('{:4d}'.format(ch))
            #    print "ch: ", '{:4d}'.format(ch),
            log_line = "ch: " + '{:4d}'.format(ch)
            for d in range(32):
                trimfile.write('{:5d}'.format(trim_final[ch][d]))
                #      print '{:5d}'.format(trim_final[ch][d]),
                log_line = log_line + '{:5d}'.format(trim_final[ch][d])
                #    print""
            log.info( log_line )
            trimfile.write("\n")

        trimfile.close()
        return 0


    def vrefpn_scan(self, pol, test_ch, npulses, amp_cal_min, amp_cal_max, amp_cal_fast, vref_t):
        # This function returns the value of Vrefp, n and Thr2_glb after a scan in a test_ch
        # at injected pulses of amplitude equal to  amp_cal_min/max and fast
        
        log.info(amp_cal_min)
        log.info(amp_cal_max)
        log.info(amp_cal_fast)
        grp_sel = test_ch%4
        
        vrefp_min = 40
        vrefp_max = 60
        ivrefp = 0
        vrefn_min = 15
        vrefn_max = 40
        ivrefn = 0
        thr2g_min = 10
        thr2g_max = 70
        ithr2g = 0
        
        diff_0 = 100
        vrefp_final = -1
        vrefn_final = -1
        thr2g_final = -1
        
        vrefp_cnt = [0 for vrefp in range(vrefp_min, vrefp_max+1)]
        vrefn_cnt = [0 for vrefn in range(vrefn_min, vrefn_max+1)]
        thr2g_cnt = [0 for thr2g in range(thr2g_min, thr2g_max+1)]

        print (" ")
        print (" ................ VRef_P,N Scan function ................")
        print (" ")
        self.write(130, 18, vref_t)
        print ("VRef_T:\t{}".format(self.read(130,18)&0xff))
        print(" ")
        print ("ASIC Efuse {} \t HW Address: {} \t Polarity {}".format(self.read_efuse_str(), self.address, pol))
        self.write(130, 5, grp_sel)
        print ("Selected channel: ", test_ch," and group ", self.read(130,5)&0xff , "\n")
        print (" ........................................................")

        if (amp_cal_min < 0. or amp_cal_min > 250.):
            amp_cal_min = 30.
            log.info("Amp_cal_min corrected to value: " + '{:3f}'.format(amp_cal_min))
        if (amp_cal_max < 0. or amp_cal_max > 250. or amp_cal_min > amp_cal_max):
            amp_cal_max = 247.
            log.info("Amp_cal_max corrected to value: " + '{:3f}'.format(amp_cal_max))

        step = (amp_cal_max - amp_cal_min)/31.
        vp_min = int(amp_cal_min)
        vp_max = int (amp_cal_max - step)

        print (" ")
        print ("Acquiring data for VRef_P scan ")
        # ..... SCANNING VRef_P ......                                                                                   
        self.write(130, 4, vp_min)
        print("vpulse set to amp_cal_min: ", self.read(130,4)&0xff)
        #self.write(130, 5, grp_shslow)
        for vrefp in range(vrefp_min, vrefp_max+1, 1):
            self.write(130, 9, vrefp & 0xff)
            self.reset_adc_counter()
            self.gen_test_pulses(n_pulses, vp_min, grp_sel)
            cnt_val = 0
            ## read counters                                                                                      
            cnt_val = self.read(test_ch,60)&0xfff
            #print(cnt_val)
            vrefp_cnt[ivrefp] = cnt_val
            ivrefp+=1
            
        print(" ")
        print("Acquiring data for VRef_N scan ")
        # ..... SCANNING VRefN .....                                                                                     
        self.write(130, 4, vp_max)
        print("vpulse set to amp_cal_max: ", self.read(130,4)&0xff)
        #self.write(130, 5, grp_shslow)
        #print "Selected channel: ", test_ch," and group ",self.read(130,5)&0xff , "\n"                              
        for vrefn in range(vrefn_min, vrefn_max+1, 1):
            self.write(130, 8, vrefn & 0xff)
            self.reset_adc_counter()
            self.gen_test_pulses(n_pulses, vp_max, grp_sel)
            cnt_val = 0
            ## read counters
            cnt_val = self.read(test_ch,0)&0xfff
            #print(cnt_val)
            vrefn_cnt[ivrefn] = cnt_val
            ivrefn+=1

        print(" ")
        print("Acquiring data for Thr2_glb scan")
        # ..... SCANNING Thr2_glb .....                                                                                  
        self.write(130, 4, vp_min)
        print("vpulse set to amp_cal_fast: ", self.read(130,4)&0xff)
        #self.write(130, 5, grp_shslow)
        #print "Selected channel: ", test_ch," and group ",self.read(130,5)&0xff , "\n"                     
        for thr2g in range(thr2g_max, thr2g_min-1, -1):
            self.write(130, 7, thr2g & 0xff)
            self.reset_adc_counter()
            self.gen_test_pulses(n_pulses, vp_min, grp_sel)
            cnt_val = 0
            ## read counters                                                                                         
            cnt_val = self.read(test_ch,62)&0xfff
            #print(cnt_val)
            thr2g_cnt[ithr2g] = cnt_val
            ithr2g+=1

        ##................... ANALYSIS ........................                                                          
        ivrefp = 0
        diff = diff_0
        for vrefp in range(vrefp_min+1, vrefp_max, 1):
            hh = int(n_pulses/2)
            if (abs(hh - vrefp_cnt[ivrefp])<diff and (vrefp_cnt[ivrefp-1] < vrefp_cnt[ivrefp]) and (vrefp_cnt[ivrefp+1]>= vrefp_cnt[ivrefp])):
                diff = abs(hh - vrefp_cnt[ivrefp])
                vrefp_final = vrefp-1
            ivrefp+=1
            
        ivrefn = 0
        diff = diff_0
        for vrefn in range(vrefn_min+1, vrefn_max, 1):
            hh = int(n_pulses/2)
            if (abs(hh - vrefn_cnt[ivrefn])<diff and (vrefn_cnt[ivrefn-1] < vrefn_cnt[ivrefn]) and (vrefn_cnt[ivrefn+1] >= vrefn_cnt[ivrefn])):
                diff = abs(hh - vrefn_cnt[ivrefn])
                vrefn_final = vrefn-1
            ivrefn+=1

        ithr2g = 0
        diff = diff_0
        for thr2g in range(thr2g_max-1, thr2g_min, -1):
            #hh = int(2*n_pulses/3)
            hh = int(n_pulses/2)
            if (abs(hh - thr2g_cnt[ithr2g])<diff and (thr2g_cnt[ithr2g -1] < thr2g_cnt[ithr2g]) and (thr2g_cnt[ithr2g+1] >= thr2g_cnt[ithr2g])):
                diff = abs(hh - thr2g_cnt[ithr2g])
                thr2g_final = thr2g+1
            ithr2g+=1

        ##................. PUBLISHING RESULTS ..................

        ivrefp = 0
        print("")
        print(" ..... Vref_P scan results .....")
        print("Test channel: ", '{:3d}'.format(test_ch))
        for vrefp in range(vrefp_min, vrefp_max+1, 1):
            if (vrefp == vrefp_final):
                print(colored("Vref_P: ", 'red',attrs=['bold']), colored( '{:2d}'.format(vrefp),'red',attrs=['bold']), end="\t")
                print(colored('{:4d}'.format(vrefp_cnt[ivrefp]),'red',attrs=['bold']), colored( "<<<---", 'red',attrs=['bold']))
            else:
                print("Vref_P: ", '{:2d}'.format(vrefp), end="\t")
                print('{:4d}'.format(vrefp_cnt[ivrefp]))
            ivrefp+=1

        ivrefn = 0
        print("")
        print(" ..... Vref_N scan results .....")
        print("Test channel: ", '{:3d}'.format(test_ch))
        for vrefn in range(vrefn_min, vrefn_max+1, 1):
            if (vrefn == vrefn_final):
                print(colored("Vref_N: ",'red',attrs=['bold']), colored('{:2d}'.format(vrefn), 'red',attrs=['bold']), end="\t")
                print(colored('{:4d}'.format(vrefn_cnt[ivrefn]), 'red',attrs=['bold']),colored("<<<---", 'red',attrs=['bold']))
            else:
                print("Vref_N: ", '{:2d}'.format(vrefn), end="\t")
                print('{:4d}'.format(vrefn_cnt[ivrefn]))
            ivrefn+=1

        ithr2g = 0
        print("")
        print(" ..... Thr2_glb scan results .....")
        print("Test channel: ", '{:3d}'.format(test_ch))
        for thr2g in range(thr2g_max, thr2g_min-1, -1):
            if (thr2g == thr2g_final):
                print(colored("Thr2_glb: ",'red',attrs=['bold']), colored('{:2d}'.format(thr2g), 'red',attrs=['bold']), end="\t")
                print(colored('{:4d}'.format(thr2g_cnt[ithr2g]), 'red',attrs=['bold']),colored("<<<---", 'red',attrs=['bold']))
            else:
                print("Thr2_glb: ", '{:2d}'.format(thr2g), end="\t")
                print('{:4d}'.format(thr2g_cnt[ithr2g]))
            ithr2g+=1

        print("")
        #thr2g_final = 25
        if (vrefp_final ==-1):
            vrefp_final = 56
            print("Vref_P value couldn't be determined, therefore, it is recommended to use: ", colored('{:2d}'.format(vrefp_final), 'red', attrs=['bold']))
        else:
            print("Vref_P value was determined and it is recommended to be set at: ", colored('{:2d}'.format(vrefp_final), 'red', attrs=['bold']))
        if (vrefn_final ==-1):
            vrefn_final = 28
            print("Vref_N value couldn't be determined, therefore, it is recommended to use: ", colored('{:2d}'.format(vrefn_final),'red', attrs=['bold']))
        else:
            print("Vref_N value was determined and it is recommended to be set at: ", colored('{:2d}'.format(vrefn_final),'red', attrs=['bold']))
        if (thr2g_final ==-1):
            thr2g_final = thr2g_min
            print("Thr2_glb value couldn't be determined, therefore, it is recommended to use: ", colored('{:2d}'.format(thr2g_final),'red', attrs=['bold']))
        else:
            print("Thr2_glb value was determined and it is recommended to be set at: ", colored('{:2d}'.format(thr2g_final),'red',attrs=['bold']))
        
        return (vrefp_final, vrefn_final, thr2g_final)
        
        
        
    def check_trim_red(self, pscan_dir, pol, asic_id_str, disc_list: int = [5, 10, 16, 24, 30, 31],vp_min: int = 0, vp_max: int = 255, vp_step: int =1, npulses: int =100):
        ''' scans over internal test pulse amplitudes to measure an S-curve
        - settings (Vref, ...) need to be done beforehand
        - measures S-curves for a subset of discriminators 
        - writes S-curve raf data to file
        args: - path for measured pscan files
              - vp_min, vp_max: min. and max. pulse amplitude
              - sn: FEB serial number ( for file name) 
        returns:
        - array with counter values for all channels, discriminators and pulse amplitudes
        '''
        disc_list = disc_list
        grps = [0,1,2,3]        
        # 99 is a dummy FEB8 SN

        # Filename elements
        if not pscan_dir.endswith('/'):
            pscan_dir += '/'
        now = datetime.now()
        date =  now.strftime("%y%m%d_%H%M")
        hw_address = self.read(192,22)
        asic_id_str = asic_id_str
        vref_p = self.read(130,9) & 0xff
        vref_n = self.read(130,8) & 0xff
        vref_t = self.read(130,18) & 0xff
        thr2_glb = self.read(130,7) & 0xff
        filename_pscan = ( "pscan_" + date
                          + "_" + str(asic_id_str)
                          + "_HW_" + str(hw_address)
                          + "_SET_" + str(vref_p) + "_" + str(vref_n) + "_" + str(vref_t) + "_" + str(thr2_glb)
                          + "_VP_" + str(vp_min) + "_" + str(vp_max) + "_"  + str(vp_step)
                          + "_NP_" + str(npulses))
        
        # ---------------  Polarity selection --------------
        #if (self.chrg_type == 1):
        if (pol == 1):
            filename_pscan = filename_pscan + "_holes.txt"
        #elif (self.chrg_type == 0):
        elif (pol == 0):
            filename_pscan = filename_pscan + "_elect.txt"
        else:
            log.error("Undefined charge type!! Exiting")
            exit()
        filename_data = pscan_dir + filename_pscan

        myfile = open(filename_data, "w+")
        myfile.write("# VP: \t CH: \t DISC_LIST:{} \t POL: {} \n".format(disc_list, pol))
        ch_step = 1
        ch_min = 0
        ch_max = 128
        ch_list = [0 for ch in range(int(ch_max/ch_step))]
        vcnt = [[[0 for vp in range(int(vp_max/vp_step))] for d in range(32)] for ch in range(int(ch_max/ch_step))]

        log.info("START ACQUIRING DATA ...")
        ivp = 0
        for vp in range(vp_min,vp_max,vp_step):
            a = (vp*15.)/256.
            # print a# Pulse height value in fC                                                                                                                        
            # print g"vp: ", '{:4d}'.format(vp), " (", '{:3f}'.format(a), ")", "\n"
            if ((vp<vp_min) or (vp>vp_max)):
                log.info(f"Pulse amplitude should be in range: {vp_min} to {vp_max}. Currently it is set to {vp}")
            ##x_data[ivp] = vp
            max_ch_grp =0
            self.reset_adc_counter()
            self.gen_test_pulses(npulses, vp, grps)
            # readback counter values
            vcnt_tmp = self.read_adc_counter(disc = disc_list, sprint = False)
            for ch in range(0,128):
                pstr = "\t"
                pstr = pstr + str(ch) +  "\t\t"
                #for d in disc_list:
                for d in disc_list:
                #for ch in range(0,128):
                    vcnt[ch][d][ivp] = vcnt_tmp[ch][d]
                    #print('{:4d}'.format(vcnt[ch][d][ivp])),
                    if d == 31:
                        pstr = pstr + "\t"
                    pstr = pstr + str(vcnt[ch][d][ivp]) + "\t"
                log.info(f"{pstr}")
            ivp+=1   

        ivp = 0
        for vp in range(vp_min,vp_max,vp_step):
            ch_counter = 0
            for ch in range(ch_min,ch_max, ch_step):
                log.info(f"vp {vp}   ch: {ch}: ")
                myfile.write("vp"+'{:4d}'.format(vp)+"   ch "+'{:4d}'.format(ch)+": ")
                #disc_list_rev = disc_list[len(disc_list)-1::-1]
                #log.info(f"{type(disc_list_rev)} {disc_list_rev}")
                #for d in disc_list_rev[1:]:      # reversed(disc_list)
                for d in disc_list:      # reversed(disc_list)
                    #print '{:4d}'.format(vcnt[ch][d_counter][ivp]),
                    myfile.write('{:6d}'.format(vcnt[ch][d][ivp]))
                ch_counter+=1
                myfile.write("\n")
            ivp+=1
        myfile.close()
        return filename_pscan

    def read_scurve_from_pscan(self, filename, disclist=[5,10,15,20,25,30]):
        # vcnt  channel disc  
        # vcnt[128][32][256]
        myfile = open(filename, "r")
        vcnt = [[[0 for vp in range(256)] for disc in range(32)] for ch in range(128)]
    
        lines = myfile.readlines()
        myfile.close()
        
        for l in lines:
            # vp 185   ch  115:      0     0    55    10   100   100
            ele = l.split()
            vp = int(ele[1])
            ch = int(ele[3][:-1])
            for d in range(4,len(ele)-1):
                disc = disclist[d-4]
                vcnt[ch][disc][vp]
        return vcnt

    def connection_check(self, pscan_dir, pol, nloops, vref_t):
        ch_min = 0
        ch_max = 128
        d_min = 0
        d_max = 32
        loop_max = nloops
        vref_t = vref_t
        # Setting the running Vref_T value
        self.write(130,18, vref_t)
        print("ASIC Efuse {}, VRef_T value {}".format(self.read_efuse_str(), self.read(130,18)&0xff))
        
        
        # Building the filename
        if not pscan_dir.endswith('/'):
            pscan_dir += '/'
        now = datetime.now()
        date =  now.strftime("%y%m%d_%H%M")
        hw_address = self.address
        asic_id_str = self.read_efuse_str()
        vref_p = self.read(130,9) & 0xff
        vref_n = self.read(130,8) & 0xff
        vref_t = self.read(130,18) & 0xff
        thr2_glb = self.read(130,7) & 0xff
        filename_conn = ( "conn_check_" + date
                          + "_" + str(asic_id_str)
                          + "_HW_" + str(hw_address)
                          + "_SET_" + str(vref_p) + "_" + str(vref_n) + "_" + str(vref_t) + "_" + str(thr2_glb)
                          + "_NL_" + str(nloops))

        if (pol == 1):
            filename_conn = filename_conn + "_holes.txt" 
        elif (pol == 0):
            filename_conn = filename_conn + "_elect.txt"
        else:
            log.error("Undefined charge type!! Exiting")
            exit()

        filename_data = pscan_dir + filename_conn
            
        myfile_conn=open(filename_data, "w+")
        vcnt = [[[0 for np in range(nloops)] for d in range(d_min,d_max)] for ch in range(ch_min,ch_max)]
        cnt_f = [[0 for d in range(d_min,d_max)] for ch in range(ch_min,ch_max)]
        # Setting the vp = 0
        #self.write(130,11, 128)
        # injecting npulses
        grp_min = 0
        grp_max = 4
        #self.gen_test_pulses(100, 80, grps)
        # Setting injected pulses amplitude to 0
        self.write(130,4,0)
        # Opening the global gate
        self.write(130,11,0)
        # running over the number of loops to read counters
        for nloops in range(0,loop_max):
            self.write(192,2,63)
            time.sleep(0.01)
            self.write(192,2,0)
            for ch in range(ch_min, ch_max):
                print('ch: {}'.format(ch), end = "\t")
                cnt_val = 0
                for d in range(d_min,d_max):
                    count = d*2
                    #print "d: ", d, " counter: ", count, "\n"                                                                     
                    cnt_val= self.read(ch,count) & 0xfff
                    cnt_f[ch][d] += cnt_val 
                    print("{:4d}".format(cnt_val), end = "\t")
                print("")

        for ch in range(ch_min,ch_max):
            myfile_conn.write("ch ")
            myfile_conn.write('{:4d}'.format(ch))
            myfile_conn.write(": ")
            for d in range(d_min,d_max):
                myfile_conn.write('{:6d}'.format(cnt_f[ch][d]))
                #print "\n"                                                                                                               
            myfile_conn.write("\n")
       
        return 0
