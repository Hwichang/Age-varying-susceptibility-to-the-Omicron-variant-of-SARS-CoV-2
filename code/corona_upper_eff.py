import numpy as np        
import pandas as pd
import copy
import random
import math
import time
###################### LOAD BASE DATA #########################################
skage_groups_new = np.array(pd.read_csv('korea_population.csv').x)
corona_new = pd.read_csv('corona_daily_omicron.csv', encoding='cp949')

##################### SCHOOL RATE ############################################
school_rate = pd.read_csv('school_rate_omi.csv')
school_rate = school_rate.T

pd.to_datetime(school_rate.columns)[0]
pd.date_range(pd.to_datetime(school_rate.columns)[0], periods=7)

school_rate_res = pd.DataFrame()

k=0
for i in range(school_rate.shape[0]-1):
    school_rate_res.loc[k,0] = pd.to_datetime(school_rate.index[i])
    school_rate_res.loc[k,1] = school_rate.iloc[i,0]
    school_rate_res.loc[k,2] = school_rate.iloc[i,1]
    school_rate_res.loc[k,3] = school_rate.iloc[i,2]
    school_rate_res.loc[k,4] = school_rate.iloc[i,3]

    k = k+1
    temp = pd.to_datetime(school_rate.index)[i]-pd.to_datetime(school_rate.index)[i+1]
    temp = temp.days
     
    if temp > 1 :
        for j in range(1,temp):
            school_rate_res.loc[k,0] = pd.to_datetime(school_rate.index)[i]-pd.Timedelta(days=j)
            school_rate_res.loc[k,1] = (1-j/temp)*(school_rate.iloc[i,0]) + (j/temp)*(school_rate.iloc[(i+1),0])
            school_rate_res.loc[k,2] = (1-j/temp)*(school_rate.iloc[i,1]) + (j/temp)*(school_rate.iloc[(i+1),1])
            school_rate_res.loc[k,3] = (1-j/temp)*(school_rate.iloc[i,2]) + (j/temp)*(school_rate.iloc[(i+1),2])
            school_rate_res.loc[k,4] = (1-j/temp)*(school_rate.iloc[i,3]) + (j/temp)*(school_rate.iloc[(i+1),3])
            k = k+1


skage = pd.read_csv('skage.csv', thousands = ',').iloc[:,4]
skage_0_2 = sum(skage[0:3])
skage_3_4 = sum(skage[3:5])
skage_5_6 = sum(skage[5:7])
skage_7_9 = sum(skage[7:10])
skage_10_12 = sum(skage[10:13])

skage_10_11 = sum(skage[10:12])
skage_12_14 = sum(skage[12:15])

skage_13_14 = sum(skage[13:15])
skage_15 = skage[15]
skage_16_18 = sum(skage[16:19])

skage_15_17 = sum(skage[15:18])
skage_18_19 = sum(skage[18:20])


skage_19 = skage[19]
skage_18_19 = sum(skage[18:20])
skage_20_29 = sum(skage[20:30])

school = np.zeros((500,4))
school_rate_res.loc[103:0:-1,1:4]
school[355:459,:] = np.array(school_rate_res.loc[103:0:-1,1:4])

################# CONTACT MATRIX ############################################
contact_school = np.array(pd.read_csv('contact_school.csv',header=None))
contact_all = np.array(pd.read_csv('contact_all.csv',header=None))
contact_work = np.array(pd.read_csv('contact_work.csv',header=None))

def contact_matrix(school_contact,S):
    res = copy.deepcopy(school_contact)
    res[0,] = (1- school[S,0]/100) * res[0,]
    res[1,] = (skage_5_6/(skage_5_6+skage_7_9))*(1-school[S,0]/100) * res[1,] +  (skage_7_9/(skage_5_6+skage_7_9))*(1-school[S,1]/100) * res[1,]
    res[2,] = (skage_10_12/(skage_10_12+skage_13_14))*(1-school[S,1]/100) * res[2,] +  (skage_13_14/(skage_10_12+skage_13_14))*(1-school[S,2]/100) * res[2,]
    res[3,] = (skage_15/(skage_15+skage_16_18+skage_19))*(1-school[S,2]/100) * res[3,] +  (skage_16_18/(skage_15+skage_16_18+skage_19))*(1-school[S,3]/100) * res[3,]
  
    res1 = contact_all - res - (0.09)*contact_work
  
    return res1


#######################InCUBATION PERIOD################
Incu_param1 = 4.544
Incu_param2 = 1/0.709


######################TRANSMISSION ONSET#################
tran_dist_mu = -4
tran_param1 = 5.2662158
tran_param2 = 1/0.8709042 


##################Infection to Recover#####################
I_R_param1 = 4
I_R_param2 = 4/5


##################Infection to Quarantine######################
C_param = 1.7


#####################Symptom to Quarantine#####################
symp_q_dist = pd.read_csv('symp_q_dist.csv').x


##################### initial SEIQ #################################
seiq_matrix = pd.read_csv('initial_seiq_omi.csv')
seiq_matrix = seiq_matrix.iloc[:,1:]


##################### initial Susceptable #################################
suscept = np.array(pd.read_csv('suscept_omi_upper_eff.csv'))
suscept_new = copy.deepcopy(suscept)


###############initial Exposed###############################
exposed =  np.array(pd.read_csv('exposed_omi.csv'))


##############initial Infectious#############################
I_total_sym =  np.array(pd.read_csv('I_total_sym_omi.csv'))
I_total_asym =  np.array(pd.read_csv('I_total_asym_omi.csv'))

I_total = I_total_sym + 0.5*I_total_asym
I_total_new = copy.deepcopy(I_total)


####################### q_matrix ###############################

# naive estimate
def naive_p(t) :
    res = (1/(skage_groups_new)) * (contact_matrix(contact_school,t-1)@I_total[(t-1),])
    
    return(res)


init_q = np.zeros((2,16))
naive_q = (exposed[1:,]/suscept[0:-1,]) / np.array(list(map(lambda x:naive_p(x), range(1,exposed.shape[0]))))


q_standard = np.array([365,439,365+100])
for i in range(1):
    temp = naive_q[q_standard[i]:(q_standard[i+1]),]
    init_q[i,:] = temp.mean(axis=0)

    
init_q[0,:] =  np.array(list(map(lambda x:max(np.random.gamma(0.001,1/0.001,1),1e-300),range(16))),dtype='float')  

q_matrix = np.zeros((500,16))

for i in range(2):
    for j in range(q_standard[i],q_standard[i+1]-1):
        q_matrix[j,] = init_q[i,:]


########################### p(t)#############################
## p_unit : force of infection at time t
## p_unit_new : force of infection at time t divided by q
## tilde_p_unit : force of infection at time t about newly imputed value
## tilde_p_unit_new : force of infection at time t divided by q newly imputed value

def p_unit(t):
    res = (q_matrix[(t-1),]/(skage_groups_new)) * (contact_matrix(contact_school,t-1)@I_total[t-1,])
    return(res)

def p_unit_new(t):
    res = (1/(skage_groups_new)) * (contact_matrix(contact_school,t-1)@I_total[t-1,])
    return(res)

def tilde_p_unit(t):
    res = (q_matrix[(t-1),]/(skage_groups_new)) * (contact_matrix(contact_school,t-1)@I_total_new[t-1,])
    return(res)

def tilde_p_unit_new(t):
    res = (1/(skage_groups_new)) * (contact_matrix(contact_school,t-1)@I_total_new[t-1,])
    return(res)  

################ For caculate P(E=t) ##########################
def p(s,t):
    res = -np.array(list(map(lambda x:p_unit(x),range(s,t+1)))).sum(axis=0)
    return(res)

def p_new(s,t):
    res = -np.array(list(map(lambda x:p_unit_new(x),range(s,t+1)))).sum(axis=0)
    return(res)

def tilde_p(s,t):
    res = -np.array(list(map(lambda x:tilde_p_unit(x),range(s,t+1)))).sum(axis=0)
    return(res)

def tilde_p_new(s,t):
    res = -np.array(list(map(lambda x:tilde_p_unit_new(x),range(s,t+1)))).sum(axis=0)
    return(res)


start_date = 366 # 2022-01-01
end_date = 397 #2022-01-31

seiq_matrix= np.array(seiq_matrix)
start = np.where(seiq_matrix[:,1]>=start_date)[0][0]
end = np.where(seiq_matrix[:,1]<=end_date+15)[0][-1]
end-start

q_list =np.zeros((1100,16))

k=1


asym_num = len(np.where(np.isnan(seiq_matrix[start:(end+1),3]))[0])
sym_num =  len(np.where(~np.isnan(seiq_matrix[start:(end+1),3]))[0])

for r in range(1100):
    sym_q = np.random.choice(np.array(symp_q_dist), size=sym_num, replace=True)
    e_sym = np.random.gamma(Incu_param1,1/Incu_param2,sym_num)
    sym_i = np.random.gamma(tran_param1,1/tran_param2,sym_num)
    asym_i = np.random.exponential(C_param,asym_num)
    asym_e = np.random.gamma(tran_param1,1/tran_param2,asym_num)
    asym_sym =np.random.gamma(Incu_param1,1/Incu_param2,asym_num)
    sym=0
    asym=0
    print(r)
    start_time = time.time()

    ######################################### - W update - ##########################################################
    for i in range(start,end+1):
        Q_i = int(seiq_matrix[i,1])-1
        i_age = int(seiq_matrix[i,0]/5)
        E_old = int(seiq_matrix[i,4])-1
        I_old = int(seiq_matrix[i,2])-1
        alpha = 0


        if ~np.isnan(seiq_matrix[i,3]) :
            c = 1
            Y_new = int(Q_i - sym_q[sym])-1
            E_new = int(Y_new -  math.ceil(e_sym[sym]))-1
            I_new = int(max(Y_new + min(math.ceil(tran_dist_mu + sym_i[sym]),14),E_new))-1
            sym += 1

        else :
            c = 0.5
            Y_new = np.nan
            I_new = int(Q_i - math.ceil(asym_i[asym]))-1
            E_new = int(I_new - max(math.ceil(asym_sym[asym] + tran_dist_mu + asym_e[asym]),0))-1
            asym += 1
        
        
        exposed[E_old,i_age] -= 1

        m_E = min(E_old,E_new)
        M_E = max(E_old,E_new)

        m_I = min(I_old,I_new)
        M_I = max(I_old,I_new)

     

        if E_old > E_new :
            suscept_new[(E_new):(E_old),i_age] -= 1

        elif E_old < E_new :
            suscept_new[(E_old):(E_new),i_age] += 1



        if I_old >= Q_i and I_new < Q_i:
            I_total_new[(I_new):(Q_i),i_age] += c  

        elif I_old > I_new and I_old < Q_i:
            I_total_new[(I_new):(I_old),i_age] += c

        elif I_old < I_new and I_new < Q_i :
            I_total_new[(I_old):(I_new),i_age] -= c

        elif I_new >= Q_i and I_old < Q_i :
            I_total_new[(I_old):(Q_i),i_age] -= c


        ##################### W_(-i) MCMC ratio ############################
        if m_I!=M_I:

            age_num_B = suscept_new[M_I,:]

            temp_B = (tilde_p(m_I+1,M_I) - p(m_I+1,M_I))*age_num_B
            alpha = alpha + sum(temp_B)
            age_num_A = exposed[(m_I+1):(M_I+1),]

            temp_A = np.array(list(map(lambda y:np.log(tilde_p_unit(y)) + tilde_p(m_I,y-1) - np.log(p_unit(y)) - p(m_I,y-1),range(m_I+1,(M_I+1)))))


            alpha = alpha + (age_num_A*temp_A).sum()
         

        ######################### W_i MCMC ratio #########################################################################


        alpha = alpha + (np.log(tilde_p_unit(E_new)) + tilde_p(m_E-1,E_new-1) - np.log(p_unit(E_old)) - p(m_E-1,E_old-1))[i_age]
        
        
        ########################### accept or reject ###################################################################

        accept = np.random.uniform(0,1,1)
        if np.log(accept) > alpha or np.isnan(alpha)==True : #reject

            I_total_new = copy.deepcopy(I_total)
            suscept_new = copy.deepcopy(suscept)
            

            exposed[E_old,i_age] += 1

        else : #accpet

            seiq_matrix[i,2] = I_new+1
            seiq_matrix[i,3] = Y_new+1
            seiq_matrix[i,4] = E_new+1

            I_total = copy.deepcopy(I_total_new)
            suscept = copy.deepcopy(suscept_new)
            exposed[E_new,i_age] += 1
    

           
    ########################################### - q update - ##########################################################
    data_num = exposed[start_date-1:end_date,]
   
    data_vec = (np.array(list(map(lambda x:-p_new(start_date-1,x) ,range(start_date-1,end_date))))*data_num).sum(axis=0)
    data_vec = data_vec - p_new(start_date-1,end_date-1)*suscept[end_date-1,]

    q_new = np.array(list(map(lambda x:np.random.gamma(0.001 + data_num.sum(axis=0)[x], 1/(0.001 + data_vec[x]),1),range(16))))
    
    q_list[k,:] = q_new[:,0]

    q_matrix[start_date-1:end_date,] = np.array(list(map(lambda x : q_new[:,0] , range(start_date-1,end_date))))

    print(q_matrix[start_date-1,:])
    pd.DataFrame(q_list).to_csv('./corona/corona_upper_eff/q_list.csv')
    pd.DataFrame(seiq_matrix[start:end+1,]).to_csv('./corona/corona_upper_eff/seiq_matrix.csv')
    pd.DataFrame(suscept).to_csv('./corona/corona_upper_eff/suscept.csv')
    pd.DataFrame(exposed).to_csv('./corona/corona_upper_eff/exposed.csv')
    pd.DataFrame(I_total).to_csv('./corona/corona_upper_eff/I_total.csv')
    k = k+1 
    
    end_time = time.time()
    print(end_time-start_time)
