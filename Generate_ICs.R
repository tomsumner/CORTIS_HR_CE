# function to generate initial conditions accounting for initial distribution of PT
#
# Takes vector initial_cond (from ODE model run), PT_S,PT_L,PT_I and E as inputs
# Returns vector of intitial conditions

IC_gen <- function(initial_cond,PT_S,PT_L,PT_I,E){

  init <- c(S=initial_cond[1]*(1-PT_S),
            S_on_P=initial_cond[1]*PT_S*E,
            S_P=initial_cond[1]*PT_S*(1-E),
            L=initial_cond[22]*(1-PT_L),
            L20=initial_cond[21]*(1-PT_L),
            L19=initial_cond[20]*(1-PT_L),
            L18=initial_cond[19]*(1-PT_L),
            L17=initial_cond[18]*(1-PT_L),
            L16=initial_cond[17]*(1-PT_L),
            L15=initial_cond[16]*(1-PT_L),
            L14=initial_cond[15]*(1-PT_L),
            L13=initial_cond[14]*(1-PT_L),
            L12=initial_cond[13]*(1-PT_L),
            L11=initial_cond[12]*(1-PT_L),
            L10=initial_cond[11]*(1-PT_L),
            L9=initial_cond[10]*(1-PT_L),
            L8=initial_cond[9]*(1-PT_L),
            L7=initial_cond[8]*(1-PT_L),
            L6=initial_cond[7]*(1-PT_L),
            L5=initial_cond[6]*(1-PT_I),
            L4=initial_cond[5]*(1-PT_I),
            L3=initial_cond[4]*(1-PT_I),
            L2=initial_cond[3]*(1-PT_I),
            L1=initial_cond[2]*(1-PT_I),
            I=0,
            on_P=E*(sum(initial_cond[7:22])*PT_L + sum(initial_cond[2:6])*PT_I),
            P=0,
            L_P=initial_cond[22]*PT_L*(1-E),
            L20_P=initial_cond[21]*PT_L*(1-E),
            L19_P=initial_cond[20]*PT_L*(1-E),
            L18_P=initial_cond[19]*PT_L*(1-E),
            L17_P=initial_cond[18]*PT_L*(1-E),
            L16_P=initial_cond[17]*PT_L*(1-E),
            L15_P=initial_cond[16]*PT_L*(1-E),
            L14_P=initial_cond[15]*PT_L*(1-E),
            L13_P=initial_cond[14]*PT_L*(1-E),
            L12_P=initial_cond[13]*PT_L*(1-E),
            L11_P=initial_cond[12]*PT_L*(1-E),
            L10_P=initial_cond[11]*PT_L*(1-E),
            L9_P=initial_cond[10]*PT_L*(1-E),
            L8_P=initial_cond[9]*PT_L*(1-E),
            L7_P=initial_cond[8]*PT_L*(1-E),
            L6_P=initial_cond[7]*PT_L*(1-E),
            L5_P=initial_cond[6]*PT_I*(1-E),
            L4_P=initial_cond[5]*PT_I*(1-E),
            L3_P=initial_cond[4]*PT_I*(1-E),
            L2_P=initial_cond[3]*PT_I*(1-E),
            L1_P=initial_cond[2]*PT_I*(1-E),
            had_TB=0)

  return(init)

}
  
  