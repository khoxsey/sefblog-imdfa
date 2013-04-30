# Copyright: Marc Wildi
# 30.07.2012
# http://blog.zhaw.ch/idp/sefblog








# New 2012 code: is needed for implementing spectral matrix in new parametrization including regularization
spec_mat_comp<-function(weight_func,L,Lag)
{
  K<-length(weight_func[,1])-1
  weight_h<-weight_func
# Frequency zero receives half weight
  weight_h[1,]<-weight_h[1,]*0.5
# Extract DFT target variable (first column)
  weight_target<-weight_h[,1]
# Rotate all DFT's such that weight_target is real (rotation does not alter mean-square error)
  weight_h<-weight_h*exp(-1.i*Arg(weight_target))
  weight_target<-weight_target*exp(-1.i*Arg(weight_target))
# DFT's explaining variables (target variable can be an explaining variable too)
  weight_h_exp<-as.matrix(weight_h[,2:(dim(weight_h)[2])])
  spec_mat<-as.vector(t(as.matrix(weight_h_exp[1,])%*%t(as.matrix(rep(1,L)))))
  for (j in 1:(K))
  {
    omegak<-j*pi/K
    exp_vec<-exp(1.i*omegak*((0:(L-1))-Lag))
    spec_mat<-cbind(spec_mat,as.vector(t(as.matrix(weight_h_exp[j+1,])%*%t(as.matrix(exp_vec)))))
  }
  dim(spec_mat)
  return(list(spec_mat=spec_mat))
}






# Modifications 30.07.2012 (new parametrization of filter coefficients)
#   1. -The asymmetry entailed by the central-deviance parametrization is avoided by proposing a new parametrization of the cross-sectional regularization
#      -The function mat_func is affected only
#      -Search for 30.07.2012 in the code below
#   2. - I removed all spectral estimation functions because they do not strictly belong to the estimation algorithm: after talking with Kent Hoxsey at JSM/San Diego we'll make a R-package where the spectral estimates (DFT,model-based,HP-based,max-entropy) will be in a separate R-file.
#   3. -the function MS_decomp_total at the end provides a decomposition of the MSE-norm into Accuracy, Timeliness and Smoothness error components which are discussed in McElroy/Wildi (2012) : you don't necessarily need this decomposition but the new function definitly belongs to this estimation-file

# Modification 06.08.2012: generalizes previous code:
#    A new parameter grand_mean is introduced
#    If grand_mean==T then previous grand-mean parametrization results (code prior to 30.07.2012)
#    If grand_mean==F then new parametrization results (from 30.07.2012 on)

# Modifications 07.08.2012 (regularization) :
#   1. We expand the effect of lambda_decay into two terms: the rate of decay (lambda_decay[1]) and the strength of the regularization (lambda_decay[2])
#   2. We disentangle L from the strength of the regularization by normlizing all matrices by their traces

# Modifications 09.08.2012 (regularization) :
# We normalize terms multiplied by des_mat in order to disentangle i1/i2 effects on the regularization troika

# 13.08.2012
# We parametrize the weights assigned to the regularization troika in such a way that the value 0 is zero_weight (no regularization) and that 1 is maximum weight (full shrinkage i.e. data is irrelevant). This is achieved with the tan-function which maps [0,pi] to the positive real numbers
# For lambda_decay[1] (shape of regularization) we impose a positive value smaller or equal than 1

# 16.08.2012: regularization works for univariate case too (in the inivariate case some matrices/vectors were not defined in the previous code)




mat_func<-function(i1,i2,L,weight_h_exp,lambda_decay,lambda_cross,lambda_smooth,Lag,weight_constraint,shift_constraint,grand_mean)
{
  if (Lag>(L-1)/2)
  {
    print("Lag larger than L/2!!!!! Will be trimmed automtically to L/2 (symmetric filter)")
    Lag<-as.integer(L/2)
  }
# 13.08.2012
# We parametrize the weights assigned to the regularization troika in such a way that the value 0 is zero_weight (no regularization) and that 1 is maximum weight (full shrinkage i.e. data is irrelevant)
# For lambda_decay[1] (shape of regularization) we impose a positive value smaller or equal than 1

  lambda_smooth<-100*tan(min(abs(lambda_smooth),0.999999)*pi/2)
  lambda_cross<-100*tan(min(abs(lambda_cross),0.999999)*pi/2)
  lambda_decay<-c(min(abs(lambda_decay[1]),1),100*tan(min(abs(lambda_decay[2]),0.999999)*pi/2))

# The smoothness and decay regularization are conveniently (rightly) implemented on original parameters
# The Q_smooth and Q_decay matrices address regularizations for original unconstrained parameters (Therefore dimension L^2)
# At the end, the matrix des_mat is used to map these regularizations to central-deviance parameters
# accounting for first order constraints!!!!!!!!!!!!!!!!!!!
  Q_smooth<-matrix(data=rep(0,((L)*length(weight_h_exp[1,]))^2),nrow=(L)*length(weight_h_exp[1,]),ncol=(L)*length(weight_h_exp[1,]))
  Q_decay<-matrix(data=rep(0,((L)*length(weight_h_exp[1,]))^2),nrow=(L)*length(weight_h_exp[1,]),ncol=(L)*length(weight_h_exp[1,]))
# Cross-sectional regularization if dimension>1
  if ((length(weight_h_exp[1,])>1))
  {
# The cross-sectional regularization is conveniently implemented on central-deviance parameters. The regularization is expressed on the
# unconstrained central-deviance parameters (dimension L), then mapped to the original (unconstrained) parameters (dimension L) with Q_centraldev_original
# and then maped back to central-deviance with constraint (dim L-1) with des_mat (mathematically unnecessarily complicate but more convenient to implement in code).
    Q_cross<-matrix(data=rep(0,((L)*length(weight_h_exp[1,]))^2),nrow=(L)*length(weight_h_exp[1,]),ncol=(L)*length(weight_h_exp[1,]))
    Q_centraldev_original<-matrix(data=rep(0,((L)*length(weight_h_exp[1,]))^2),nrow=(L)*length(weight_h_exp[1,]),ncol=(L)*length(weight_h_exp[1,]))
  } else
  {
# 16.08.2012
    Q_cross<-NULL
  }
  for (i in 1:L)
  {
# For symmetric filters or any historical filter with Lag>0 the decay must be symmetric about b_max(0,Lag) zunehmen

# New 07.08.2012: lambda_decay is now a 2-dim vector: the first component controls for the exponential decay and the second accounts for the strength of the regularization
    Q_decay[i,i]<-(1+lambda_decay[1])^(2*abs(i-1-max(0,Lag)))
    if(i==1)
    {
      Q_smooth[i,i:(i+2)]<-c(1,-2,1)
    } else
    {
      if(i==2)
      {
        Q_smooth[i,(i-1):(i+2)]<-c(-2,5,-4,1)
      } else
      {
        if(i==L)
        {
          Q_smooth[i,(i-2):i]<-c(1,-2,1)
        } else
        {
          if(i==L-1)
          {
            Q_smooth[i,(i-2):(i+1)]<-c(1,-4,5,-2)
          } else
          {
            Q_smooth[i,(i-2):(i+2)]<-c(1,-4,6,-4,1)
          }
        }
      }
    }
  }

  if (length(weight_h_exp[1,])>1)
  {

    for (j in 1:max(1,(length(weight_h_exp[1,])-1)))   #j<-1
    {
      Q_smooth[j*L+1:L,j*L+1:L]<-Q_smooth[1:L,1:L]
      Q_decay[j*L+1:L,j*L+1:L]<-Q_decay[1:L,1:L]
    }
    Q_centraldev_original<-diag(rep(1,L*length(weight_h_exp[1,])))
    diag(Q_centraldev_original[1:L,L+1:L])<-rep(-1,L)
    for (i in 2:length(weight_h_exp[1,]))   #i<-2
    {
      diag(Q_centraldev_original[(i-1)*L+1:L,1:L])<-rep(1,L)
      diag(Q_centraldev_original[(i-1)*L+1:L,(i-1)*L+1:L])<-rep(1,L)
      diag(Q_centraldev_original[1:L,(i-1)*L+1:L])<-rep(-1,L)
    }
    Q_centraldev_original<-solve(Q_centraldev_original)

# 06.08.2012: the following if allows for either grand-mean parametrization (I-MDFA version prior to 30.07.2012) or original parameters
#   If grand_mean==T then the code replicates I-MDFA as released prior 30.07.2012
#   If grand_mean==F then the new parametrization is used.
#   Differences between both approaches: see section 7.2 of my elements paper posted on SEFBlog (both codes are identical when no regularization is imposed. Otherwise the later version (grand_mean==F) is logically more consistent becuase it treats all series identically (no asymmetry)).
    if (grand_mean)
    {
      diag(Q_cross[L+1:((length(weight_h_exp[1,])-1)*L),L+1:((length(weight_h_exp[1,])-1)*L)])<-
      rep(1,((length(weight_h_exp[1,])-1)*L))
    } else
    {
#30.07.2012:new definition (parametrization) of Q_cross (Lambda_{cross} in the elements-paper)
      diag(Q_cross)<-1
      for (i in 1:length(weight_h_exp[1,]))
      {
        for (j in 1:L)
        {
          Q_cross[(i-1)*L+j,j+(0:(length(weight_h_exp[1,])-1))*L]<-Q_cross[(i-1)*L+j,j+(0:(length(weight_h_exp[1,])-1))*L]-1/length(weight_h_exp[1,])
        }
      }
    }
  } else
  {
# 16.08.2012: define matrix for univariate case
    Q_centraldev_original<-NULL
  }
# New 07.08.2012 : the next lines for normalizing the troika are new: disentangle the effect by L
  Q_decay<-Q_decay*lambda_decay[2]
  Q_cross<-Q_cross*lambda_cross                   #Qh<-Q_cross
  Q_smooth<-Q_smooth*lambda_smooth
  if (lambda_decay[2]>0)
  {
# The second parameter in lambda_decay accounts for the strength of the regularization
    Q_decay<-lambda_decay[2]*(Q_decay/(sum(diag(Q_decay))))
  }
  if (lambda_cross>0)
  {
    Q_cross<-lambda_cross*(Q_cross/(sum(diag(Q_cross))))
  }
  if (lambda_smooth>0)
  {
    Q_smooth<-lambda_smooth*(Q_smooth/(sum(diag(Q_smooth))))
  }

# new 21.06.2012: w_eight vector
# weight vector (relevant if i1<-T)     Lag<--2
  if (i1)
  {
    if (i2)
    {
# Modifications 17.04.2012 : the definition of the vector w_eight has been generalized:
#               The new definition allows to impose constraints to b_Lag, b_{Lag+1} instead of b_{L-1} and b_L
#               Therefore the decay regularization does not potentially conflict with filter constraints
      if (Lag<1)
      {
        w_eight<-c(-(Lag-1)*weight_constraint[1]-shift_constraint[1],
        Lag*weight_constraint[1]+shift_constraint[1],rep(0,L-2))
      } else
      {
        w_eight<-c(rep(0,Lag),weight_constraint[1]-shift_constraint[1],
        shift_constraint[1],rep(0,L-Lag-2))
      }

      if (length(weight_h_exp[1,])>1)
      {
        for (j in 2:length(weight_h_exp[1,]))
        {

          if (Lag<1)
          {
            w_eight<-c(w_eight,-(Lag-1)*weight_constraint[j]-shift_constraint[j],Lag*weight_constraint[j]+shift_constraint[j],rep(0,L-2))
          } else
          {
            w_eight<-c(w_eight,c(rep(0,Lag),weight_constraint[j]-shift_constraint[j],
            shift_constraint[j],rep(0,L-Lag-2)))
          }
        }
      }
    } else
    {
      if (Lag<1)
      {
        w_eight<-c(weight_constraint[1],rep(0,L-1))
      } else
      {
        w_eight<-c(rep(0,Lag),weight_constraint[1],rep(0,L-Lag-1))
      }
      if (length(weight_h_exp[1,])>1)
      {
        for (j in 2:length(weight_h_exp[1,]))
        {
          if (Lag<1)
          {
            w_eight<-c(w_eight,weight_constraint[j],rep(0,L-1))
          } else
          {
            w_eight<-c(w_eight,rep(0,Lag),weight_constraint[j],rep(0,L-Lag-1))
          }
        }
      }
    }
  } else
  {
    if (i2)
    {
      if (Lag<1)
      {
        w_eight<-c(0,shift_constraint[1]/(1-Lag),rep(0,L-2))
      } else
      {
        w_eight<-c(rep(0,Lag+1),shift_constraint[1],rep(0,L-Lag-2))
      }

      if (length(weight_h_exp[1,])>1)
      {
        for (j in 2:length(weight_h_exp[1,]))
        {
          if (Lag<1)
          {
            w_eight<-c(w_eight,c(0,shift_constraint[j]/(1-Lag),rep(0,L-2)))
          } else
          {
            w_eight<-c(w_eight,c(rep(0,Lag+1),shift_constraint[j],rep(0,L-Lag-2)))
          }
        }
      }
    } else
    {
      w_eight<-rep(0,L*length(weight_h_exp[1,]))
    }
  }

# Here we implement the matrix which links freely determined central-deviance parameters and constrained original parameters
# In my `elements'-paper t(des_mat) corresponds to A%*%R (R links constrained and unconstrained parameters and A maps central-deviance to original parameters)
#   Please note that:
#   1. Here I'm working with central-deviance parameters (in the paper I'm working with original parameters)
#   2. The same matrix R applies to either parameter set
#   3. If I work with central-deviance parameters then R maps the freely determined set to the constrained (central-deviance)
#       and A then maps the constrained (central-deviance) set to original constrained parameters.
# Modifications 17.04.2012: we generalize the definition of des_mat in the case of i2<-T (i2<-F is OK)
  if (i2)
  {
      if (i1)
      {
  # First and second order restrictions
        des_mat<-matrix(data=rep(0,(L-2)*L*(length(weight_h_exp[1,]))^2),nrow=(L-2)*length(weight_h_exp[1,]),ncol=(L)*length(weight_h_exp[1,]))

        for (i in 1:(L-2))
        {
          if (Lag<1)                                  #des_mat[1:23,]
          {
            des_mat[i,i+2+(0:(length(weight_h_exp[1,])-1))*L]<-1
            des_mat[i,1+(0:(length(weight_h_exp[1,])-1))*L]<-i
            des_mat[i,2+(0:(length(weight_h_exp[1,])-1))*L]<--(i+1)

          } else
          {
            des_mat[i,ifelse(i<Lag+1,i,i+2)+(0:(length(weight_h_exp[1,])-1))*L]<-1
            des_mat[i,Lag+1+(0:(length(weight_h_exp[1,])-1))*L]<-ifelse(i<Lag+1,-(Lag+2-i),i-Lag)
            des_mat[i,Lag+2+(0:(length(weight_h_exp[1,])-1))*L]<-ifelse(i<Lag+1,(Lag+1-i),-(i-Lag+1))
          }
        }
        if (length(weight_h_exp[1,])>1)
        {
          for (j in 1:max(1,(length(weight_h_exp[1,])-1)))
          {
            for (i in 1:(L-2))                                                          #reg_mat[600,600
            {

              if (Lag<1)
              {
                des_mat[i+j*(L-2),i+2]<--1
                des_mat[i+j*(L-2),1]<--i
                des_mat[i+j*(L-2),2]<-(i+1)
              } else
              {
                des_mat[i+j*(L-2),ifelse(i<Lag+1,i,i+2)]<--1
                des_mat[i+j*(L-2),Lag+1]<--ifelse(i<Lag+1,-(Lag+2-i),i-Lag)
                des_mat[i+j*(L-2),Lag+2]<--ifelse(i<Lag+1,(Lag+1-i),-(i-Lag+1))
              }

              if (Lag<1)
              {
                des_mat[i+j*(L-2),i+2+j*L]<-1
                des_mat[i+j*(L-2),1+j*L]<-i
                des_mat[i+j*(L-2),2+j*L]<--(i+1)
              } else
              {
                des_mat[i+j*(L-2),ifelse(i<Lag+1,i+j*L,i+2+j*L)]<-1
                des_mat[i+j*(L-2),Lag+1+j*L]<-ifelse(i<Lag+1,-(Lag+2-i),i-Lag)
                des_mat[i+j*(L-2),Lag+2+j*L]<-ifelse(i<Lag+1,(Lag+1-i),-(i-Lag+1))
              }
            }
          }
        }
      } else
      {
# Modifications 17.04.2012: we generalize the definition of des_mat in the case of i2<-T (i2<-F is OK)
# new 21.06.2012: new definition of des_mat in the case i1<-F and i2<-T
        des_mat<-matrix(data=rep(0,(L-1)*L*(length(weight_h_exp[1,]))^2),nrow=(L-1)*length(weight_h_exp[1,]),ncol=(L)*length(weight_h_exp[1,]))
        for (i in 1:(L-1))      #Lag<--2
        {

          if (Lag<1)         #des_mat[1:23,1:24]    #des_mat[1:23,]
          {
# For Lag<1 we have to impose the restriction to b2 (b1 does not appear in the i2<-T and i1<-F case!)
            des_mat[i,ifelse(i<2,i,i+1)+(0:(length(weight_h_exp[1,])-1))*L]<-1
#            des_mat[i,1+1+(0:(length(weight_h_exp[1,])-1))*L]<-ifelse(i<2,-(i-1),-i)/(1)
            des_mat[i,1+1+(0:(length(weight_h_exp[1,])-1))*L]<-ifelse(i==1,Lag,Lag-i)/(1-Lag)

          } else
          {
            des_mat[i,ifelse(i<Lag+2,i,i+1)+(0:(length(weight_h_exp[1,])-1))*L]<-1
            des_mat[i,Lag+2+(0:(length(weight_h_exp[1,])-1))*L]<-ifelse(i<Lag+2,Lag+1-i,Lag-i)
          }
        }
# new 08.02.2012: the if loop is new/modified
        if (length(weight_h_exp[1,])>1)
        {
          for (j in 1:max(1,(length(weight_h_exp[1,])-1)))
          {
            for (i in 1:(L-1))
            {

              if (Lag<1)
              {
    # For Lag<1 we have to impose the restriction to b2 (b1 does not appear in the i2<-T and i1<-F case!)
                des_mat[i+j*(L-1),ifelse(i<2,i,i+1)]<--1
#                des_mat[i+j*(L-1),1+1]<--ifelse(i<1+1,-(i-1),-i)/(1)
                des_mat[i+j*(L-1),1+1]<--ifelse(i==1,Lag,Lag-i)/(1-Lag)
              } else
              {
                des_mat[i+j*(L-1),ifelse(i<Lag+2,i,i+1)]<--1
                des_mat[i+j*(L-1),Lag+2]<--ifelse(i<Lag+2,Lag+1-i,Lag-i)
              }

              if (Lag<1)
              {
    # For Lag<1 we have to impose the restriction to b2 (b1 does not appear in the i2<-T and i1<-F case!)
                des_mat[i+j*(L-1),ifelse(i<2,i,i+1)+j*L]<-1
#                des_mat[i+j*(L-1),1+1+j*L]<-ifelse(i<1+1,-(i-1),-i)/(1)
                des_mat[i+j*(L-1),1+1+j*L]<-ifelse(i==1,Lag,Lag-i)/(1-Lag)
              } else
              {
                des_mat[i+j*(L-1),ifelse(i<Lag+2,i,i+1)+j*L]<-1
                des_mat[i+j*(L-1),Lag+2+j*L]<-ifelse(i<Lag+2,Lag+1-i,Lag-i)
              }

            }
          }
        }

      }
  } else
  {
    if (i1)                            #lambda_cross=0
    {
      des_mat<-matrix(data=rep(0,(L-1)*L*(length(weight_h_exp[1,]))^2),nrow=(L-1)*length(weight_h_exp[1,]),ncol=(L)*length(weight_h_exp[1,]))
# The cross-sectional regularization can be directly implemented into reg_mat (it adresses the central deviance parameters)
      for (i in 1:(L-1))  #i<-1
      {
# The i1-constraint is imposed on b_max(0,Lag) (instead of b_L as in original procedure) in order to avoid a conflict with the exponential decay requirement
        if (Lag<1)
        {
          des_mat[i,i+1+(0:(length(weight_h_exp[1,])-1))*L]<-1
          des_mat[i,1+(0:(length(weight_h_exp[1,])-1))*L]<--1
        } else
        {
# Lag cannot be larger than (L-1)/2 (symmetric filter)
          des_mat[i,ifelse(i<Lag+1,i,i+1)+(0:(length(weight_h_exp[1,])-1))*L]<-1
          des_mat[i,Lag+1+(0:(length(weight_h_exp[1,])-1))*L]<--1
        }
      }
# new 08.02.2012: the if loop is new/modified
      if (length(weight_h_exp[1,])>1)
      {
        for (j in 1:max(1,(length(weight_h_exp[1,])-1)))   #j<-1
        {
          for (i in 1:(L-1))
          {
  # The i1-constraint is imposed on b_max(0,Lag) (instead of b_L as in original procedure) in order to avoid a conflict with the exponential decay requirement
            if (Lag<1)
            {
              des_mat[i+j*(L-1),i+1]<--1
              des_mat[i+j*(L-1),1]<-1
              des_mat[i+j*(L-1),i+1+j*L]<-1
              des_mat[i+j*(L-1),1+j*L]<--1
            } else
            {
  # Lag cannot be larger than (L-1)/2 (symmetric filter)
              des_mat[i+j*(L-1),ifelse(i<Lag+1,i,i+1)]<--1
              des_mat[i+j*(L-1),Lag+1]<-1
              des_mat[i+j*(L-1),ifelse(i<Lag+1,i,i+1)+j*L]<-1
              des_mat[i+j*(L-1),Lag+1+j*L]<--1
            }
  # The cross sectional regularization is implemented directly into reg_mat (it addresses central deviance parameters!)
#            reg_mat[i+(j)*(L-1),i+(j)*(L-1)]<-lambda_cross
          }
        }             # det(reg_mat)         lambda_cross<-0.3
      }
    } else
    {
      des_mat<-matrix(data=rep(0,(L)*L*(length(weight_h_exp[1,]))^2),nrow=(L)*length(weight_h_exp[1,]),ncol=(L)*length(weight_h_exp[1,]))
      for (i in 1:(L))
      {
        des_mat[i,i+(0:(length(weight_h_exp[1,])-1))*L]<-1
      }
# new 08.02.2012: the if loop is new/modified
      if (length(weight_h_exp[1,])>1)
      {

        for (j in 1:max(1,(length(weight_h_exp[1,])-1)))
        {
          for (i in 1:(L))
          {
            des_mat[i+(j)*(L),i]<--1
            des_mat[i+(j)*(L),i+j*L]<-1
          }
        }
      }
    }
  }

# 6.08.2012: the following if allows for either grand-mean parametrization (I-MDFA version prior to 30.07.2012) or original parameters
# 16.08.2012: the additional condition length(weight_h_exp[1,])>1) is checked
  if (!grand_mean&length(weight_h_exp[1,])>1)
  {
# 30.07.2012: patch which does apply A^{-1} to A*R giving R where A and R are defined in my elememts paper
#and des_mat=A*R i.e.  t(Q_centraldev_original%*%t(des_mat)) is R (called des_mat in my code)
    des_mat<-t(Q_centraldev_original%*%t(des_mat))
  }
# Here we fold all three regularizations (cross, smooth and decay) into a single reg-matrix
# The smoothness and decay terms address original parameters and must be transformed (by des_mat) in order to
# conform to the central-deviance parameterization as well as to first (and/or second) order constraints

## begin new 08.02.2012: the if-loop is new

  if ((length(weight_h_exp[1,])>1))
  {
# 06.08.2012: the following if allows for either grand-mean parametrization (I-MDFA version prior to 30.07.2012) or original parameters
    if (grand_mean)
    {
      reg_t<-(Q_smooth+Q_decay+t(Q_centraldev_original)%*%Q_cross%*%Q_centraldev_original)
    } else
    {
#    reg_t<-(Q_smooth+Q_decay+t(Q_centraldev_original)%*%Q_cross%*%Q_centraldev_original)
# 30.07.2012 The cross-sectional regularization term is now simpler (without A-transformation)
      reg_t<-(Q_smooth+Q_decay+Q_cross)
    }
  } else
  {
    reg_t<-(Q_smooth+Q_decay)
  }

# 09.08.2012  : normalize regularization terms (which are multiplied by des_mat) in order to disentangle i1/i2 effects
  reg_mat<-(des_mat)%*%reg_t%*%t(des_mat)
  if (lambda_smooth+lambda_decay[2]+lambda_cross>0)
  {
    disentangle_des_mat_effect<-sum(diag(reg_t))/sum(diag(reg_mat))#sum(apply(reg_mat,1,sum))/sum(apply(reg_t,1,sum))
    reg_mat<-reg_mat*disentangle_des_mat_effect
    reg_xtxy<-des_mat%*%reg_t%*%w_eight*(disentangle_des_mat_effect)#+t(w_eight)%*%reg_t%*%t(des_mat)
  } else
  {
    reg_xtxy<-des_mat%*%reg_t%*%w_eight#+t(w_eight)%*%reg_t%*%t(des_mat)
  }


## end new 08.02.2012
  return(list(des_mat=des_mat,reg_mat=reg_mat,reg_xtxy=reg_xtxy,w_eight=w_eight))
}


#Q_centraldev_original%*%t(des_mat)
#solve(Q_centraldev_original)




















# New 2012 I-MDFA code: might crash if dimension is 1 (univariate setting): I did not test this!!!
# Generalizes old code (below) when regularization parameters are set to zero
# This new function can deal with very richly parametrized designs (high-dimensional multivariate with long lag structures)
# Regularization affect/control smoothness, rate of decay and cross-sectional similarity of filter parameters/weights

# new 21.06.2012: new vector shift_constraint in function call

# 06.08.2012: new parameter grand_mean in function call
mdfa_analytic_new<-function(K,L,lambda,weight_func,Lag,Gamma,expweight,cutoff,i1,i2,weight_constraint,lambda_cross,lambda_decay,lambda_smooth,lin_expweight,shift_constraint,grand_mean)
{

# In order to enhance numerical speed this call could be done outside (as long as L and Lag are fixed)
  spec_mat<-spec_mat_comp(weight_func,L,Lag)$spec_mat     #dim(spec_mat[,1])
# weighting of amplitude function in stopband
  omega_Gamma<-as.integer(cutoff*K/pi)
  if ((K-omega_Gamma+1)>0)
  {
#    lin_expweight <- FALSE
    if (lin_expweight)
    {
      expweight_vec<-c(rep(1,omega_Gamma),1+rep(expweight,K-omega_Gamma+1))
    } else
    {
      expweight_vec<-c(rep(1,omega_Gamma),(1:(K-omega_Gamma+1))^(expweight/2))
#      expweight_vec<-c(rep(1,omega_Gamma),1+(1:(K-omega_Gamma+1))*(expweight/2))
    }
    weight_h<-weight_func*expweight_vec
  } else
  {
    expweight_vec<-rep(1,K+1)
    weight_h<-weight_func* expweight_vec
  }
#ts.plot(abs(weight_h))       ts.plot(Gamma)
# Frequency zero receives half weight
  weight_h[1,]<-weight_h[1,]*0.5
# DFT target variable
  weight_target<-weight_h[,1]
# Rotate all DFT's such that weight_target is real (rotation does not alter mean-square error)
  weight_h<-weight_h*exp(-1.i*Arg(weight_target))
  weight_target<-Re(weight_target*exp(-1.i*Arg(weight_target)))
# DFT's explaining variables: target variable can be an explaining variable too
  weight_h_exp<-as.matrix(weight_h[,2:(dim(weight_h)[2])])

# The spectral matrix is inflated in stopband: effect of expweight
  spec_mat<-t(t(spec_mat)*expweight_vec) #dim(spec_mat)

# Compute design matrix and regularization matrix
# new 21.06.2012: new vector shift_constraint in function call
# 06.08.2012: new parameter grand_mean in function call
  mat_obj<-mat_func(i1,i2,L,weight_h_exp,lambda_decay,lambda_cross,lambda_smooth,Lag,weight_constraint,shift_constraint,grand_mean)

  des_mat<-mat_obj$des_mat
  reg_mat<-mat_obj$reg_mat           #dim(des_mat[1:23,1:24])
  reg_xtxy<-mat_obj$reg_xtxy
  w_eight<-mat_obj$w_eight

# Solve estimation problem in space of b_main and b_dev
  mat_x<-des_mat%*%spec_mat          #dim(mat_x)    dim(spec_mat)   length(w_eight)
  X_new<-t(Re(mat_x))+sqrt(1+Gamma*lambda)*1.i*t(Im(mat_x))
# xtx can be written either in Re(t(Conj(X_new))%*%X_new) or as below:
  xtx<-t(Re(X_new))%*%Re(X_new)+t(Im(X_new))%*%Im(X_new)
# The filter restrictions (i1<-T and/or i2<-T) appear as constants on the right hand-side of the equation:
  xtxy<-t(Re(t(w_eight)%*%spec_mat)%*%Re(t(spec_mat)%*%t(des_mat))+
  Im(t(w_eight)%*%t(t(spec_mat)*sqrt(1+Gamma*lambda)))%*%Im(t(t(t(spec_mat)*sqrt(1+Gamma*lambda)))%*%t(des_mat)))
# scaler makes scales of regularization and unconstrained optimization `similar'
  scaler<-mean(diag(xtx))
  X_inv<-solve(xtx+scaler*reg_mat)
  bh<-as.vector(X_inv%*%(((t(Re(X_new)*weight_target))%*%Gamma)-xtxy-scaler*reg_xtxy))
# the last two filter weights are functions of the previous ones through the first and second order restrictions
  b<-matrix(nrow=L,ncol=length(weight_h_exp[1,]))
# Reconstruct original parameters
  bhh<-t(des_mat)%*%bh
  for (k in 1:L) #k<-1
  {      #dim(t(des_mat))
    b[k,]<-bhh[(k)+(0:(length(weight_h_exp[1,])-1))*L]
  }

# Modifications 17.04.2012 : the newly defined vector w_eight allows for simple/straightforward adjustment of filter coefficients
  weight_cm<-matrix(w_eight,ncol=(length(weight_h_exp[1,])))
# Add level and/or time-shift constraints (if i1<-F and i2<-F then this matrix is zero)
  b<-b+weight_cm

# The following derivations of the DFA-criterion are all equivalent
# They are identical with rever (up to normalization by (2*(K+1)^2))  as calculated at the end of the function except if i1<-T and weight_constraint different from zero
# In the latter case an additional constant interfers with the estimation

# The target Y in the frequency-domain is the real vector weight_target*Gamma (both vectors are real)
# The Regression estimate (the smoother) of Y is the following expression:
  trth<-((X_new)%*%(X_inv%*%t(Re(X_new))))%*%(weight_target*Gamma)                    #sum(abs(trth-trt))
# This expression is identical to trt computed below if lambda=0 (assuming i1<-F or weight_constraint=0); otherwise trth is identical to Re(trt)+1.i*sqrt(1+lambda*Gamma)*Im(trt))
# The projection matrix (it is a projection for the real part only, see below) is therefore:
  Proj_mat<-((X_new)%*%(X_inv%*%t(Re(X_new))))              #trth-(Proj_mat)%*%trth
# The residual projection matrix is (it is a projection for the real part, see below)
  res_mat<-diag(rep(1,dim(Proj_mat)[1]))-Proj_mat
# DFA criterion: first possibility (all three variants are identical)
  sum(abs(res_mat%*%(weight_target*Gamma))^2)
# Residuals (DFT of filter error): in contrast to OLS this is not iid (weighted regression in frequency-domain)
  resi<-res_mat%*%(weight_target*Gamma)
# DFA criterion: second possibility
  t(Conj(resi))%*%resi
  t((weight_target*Gamma))%*%(t(Conj(res_mat))%*%(res_mat))%*%(weight_target*Gamma)
# The interesting `effective degrees of freedom' used here emphasizes an unbiased estimate of the mean-square residual
#    See  http://en.wikipedia.org/wiki/Degrees_of_freedom_(statistics) (effective degrees of freedom: the expression tr((I-H)(I-H)')
#    Note that res_mat=I-H see above
#    Then (Y-Y^)(Y-Y^)/tr((I-H)(I-H)') is an unbiased estimate of the mean-square residual error (in our case Y=weight_target*Gamma. see above)
#    This correcting `effective degrees of freedom' can then be used to implement a generalized AIC, see below
  degrees_freedom<-2*Re(sum(diag(t(Conj(res_mat))%*%(res_mat))))-1
# Note that res_mat is a projection matrix with regards to the real part (but not with respect to the imaginary part)
# Thus we could replace diag(t(Conj(res_mat))%*%(res_mat)) by diag(res_mat) in the degrees_freedom above i.e. the following expression equates to zero:
  Re(t(Conj(res_mat))%*%(res_mat))-Re(res_mat)
  freezed_degrees<-2*K+1-degrees_freedom
# This is an alternative (identical) expression for the freezed_degreees
  2*Re(sum(diag(Proj_mat)))

# DFA Criterion: third possibility (here an additional normalization by 1/(2*(K+1)^2))
  sum(abs(Gamma*weight_target-trth)^2)

#ts.plot(b)
  # Transferfunction
  trffkt<-matrix(nrow=K+1,ncol=length(weight_h_exp[1,]))
  trffkth<-trffkt
  trffkt[1,]<-apply(b,2,sum)
  trffkth[1,]<-trffkt[1,]
#  b<-scale(b,center=T,scale=F)

  for (j in 1:length(weight_h_exp[1,]))
  {
    for (k in 0:(K))#k<-1
    {
      trffkt[k+1,j]<-(b[,j]%*%exp(1.i*k*(0:(L-1))*pi/(K)))
    }
  }
  trt<-apply(((trffkt)*exp(1.i*(0-Lag)*pi*(0:(K))/K))*weight_h_exp,1,sum)
# DFA criterion which accounts for customization but not for regularization term
  rever<-sum(abs(Gamma*weight_target-Re(trt)-1.i*sqrt(1+lambda*Gamma)*Im(trt))^2)/(2*(K+1)^2)
# MS-filter error : DFA-criterion without effects by lambda or expweight (one must divide spectrum by expweight_vec)
  MS_error<-sum((abs(Gamma*weight_target-trt)/expweight_vec)^2)/(2*(K+1)^2)
# Definition of Accuracy, time-shift and noise suppression terms
# Please note that:
#       1. these terms refer to the original non-linearized criterion: therefore they do not sum up to rever
#       2. we are interested in decomposing the mean-square error i.e. we ignore expweight and lambda here (we just want to measure the impact of lambda and expweight)
#               Therefore, we use the DFT without expweight_vec i.e. we have to divide all spectral constents by expweight_vec
  Gamma_cp<-Gamma[1+0:as.integer(K*(cutoff/pi))]
  Gamma_cn<-Gamma[(2+as.integer(K*(cutoff/pi))):(K+1)]
  trt_cp<-(trt/expweight_vec)[1+0:as.integer(K*(cutoff/pi))]
  trt_cn<-(trt/expweight_vec)[(2+as.integer(K*(cutoff/pi))):(K+1)]
  weight_target_cp<-(weight_target/expweight_vec)[1+0:as.integer(K*(cutoff/pi))]
  weight_target_cn<-(weight_target/expweight_vec)[(2+as.integer(K*(cutoff/pi))):(K+1)]
# define singular observations
  Accuracy<-sum(abs(Gamma_cp*weight_target_cp-abs(trt_cp))^2)/(2*(K+1)^2)
  Timeliness<-4*sum(abs(Gamma_cp)*abs(trt_cp)*sin(Arg(trt_cp)/2)^2*weight_target_cp)/(2*(K+1)^2)
  Smoothness<-sum(abs(Gamma_cn*weight_target_cn-abs(trt_cn))^2)/(2*(K+1)^2)
  Shift_stopband<-4*sum(abs(Gamma_cn)*abs(trt_cn)*sin(Arg(trt_cn)/2)^2*weight_target_cn)/(2*(K+1)^2)
# Check: the following expression should vanish
  Accuracy+Timeliness+Smoothness+Shift_stopband-MS_error
# Very prototypical: AIC
#  aic<-ifelse(degrees_freedom<K+1&degrees_freedom>1,log(rever)+2*(K-degrees_freedom)/(K)+2*(K-degrees_freedom)*(K-degrees_freedom+1)/(K*(degrees_freedom-1)),NA)
  aic<-ifelse(degrees_freedom<K+1&degrees_freedom>1,log(rever)+2*(K-degrees_freedom+1)/(degrees_freedom-2),NA)

  return(list(b=b,trffkt=trffkt,rever=rever,degrees_freedom=degrees_freedom,aic=aic,freezed_degrees=freezed_degrees,Accuracy=Accuracy,Smoothness=Smoothness,Timeliness=Timeliness,MS_error=MS_error))
}












MS_decomp_total<-function(Gamma,trffkt,weight_func,cutoff,Lag)     #trffkt<-trffkt_hp         pi/ub
{
  if (!(length(trffkt[,1])==length(weight_func[,1])))
  {
    len_w<-min(length(trffkt[,1]),length(weight_func[,1]))
    if (length(trffkt[,1])<length(weight_func[,1]))
    {
      len_r<-(length(weight_func[,1])-1)/(length(trffkt[,1])-1)
      weight_funch<-weight_func[c(1,(1:(len_w-1))*len_r),]
      trffkth<-trffkt
    } else
    {
      len_r<-1/((length(weight_func[,1])-1)/(length(trffkt[,1])-1))
      trffkth<-trffkt[c(1,(1:(len_w-1))*len_r),]
      weight_funch<-weight_func
    }
  } else
  {
    len_w<-length(trffkt[,1])
    weight_funch<-weight_func
    trffkth<-trffkt
    Gammah<-Gamma
  }
  if (length(Gamma)>len_w)
  {
    len_r<-(length(Gamma)-1)/(len_w-1)
    Gammah<-Gamma[c(1,(1:(len_w-1))*len_r)]
  }

#cbind(Ghh,Gammah)

  weight_h<-weight_funch
  K<-length(weight_funch[,1])-1
  weight_target<-weight_h[,1]
# Rotate all DFT's such that weight_target is real (rotation does not alter mean-square error)
  weight_h<-weight_h*exp(-1.i*Arg(weight_target))
  weight_target<-Re(weight_target*exp(-1.i*Arg(weight_target)))
# DFT's explaining variables: target variable can be an explaining variable too
  weight_h_exp<-as.matrix(weight_h[,2:(dim(weight_h)[2])])

  trt<-apply(((trffkth)*exp(1.i*(0-Lag)*pi*(0:(K))/K))*weight_h_exp,1,sum)
# MS-filter error : DFA-criterion without effects by lambda or expweight (one must divide spectrum by expweight_vec)
  MS_error<-sum((abs(Gammah*weight_target-trt))^2)/(2*(K+1)^2)
# Definition of Accuracy, time-shift and noise suppression terms
# Please note that:
#       1. these terms refer to the original non-linearized criterion: therefore they do not sum up to rever
#       2. we are interested in decomposing the mean-square error i.e. we ignore expweight and lambda here (we just want to measure the impact of lambda and expweight)
#               Therefore, we use the DFT without expweight_vec i.e. we have to divide all spectral constents by expweight_vec
  Gamma_cp<-Gammah[1+0:as.integer(K*(cutoff/pi))]
  Gamma_cn<-Gammah[(2+as.integer(K*(cutoff/pi))):(K+1)]
  trt_cp<-trt[1+0:as.integer(K*(cutoff/pi))]
  trt_cn<-trt[(2+as.integer(K*(cutoff/pi))):(K+1)]
  weight_target_cp<-weight_target[1+0:as.integer(K*(cutoff/pi))]
  weight_target_cn<-weight_target[(2+as.integer(K*(cutoff/pi))):(K+1)]
# define singular observations
  Accuracy<-sum(abs(Gamma_cp*weight_target_cp-abs(trt_cp))^2)/(2*(K+1)^2)
  Timeliness<-4*sum(abs(Gamma_cp)*abs(trt_cp)*sin(Arg(trt_cp)/2)^2*weight_target_cp)/(2*(K+1)^2)
  Smoothness<-sum(abs(Gamma_cn*weight_target_cn-abs(trt_cn))^2)/(2*(K+1)^2)
  Shift_stopband<-4*sum(abs(Gamma_cn)*abs(trt_cn)*sin(Arg(trt_cn)/2)^2*weight_target_cn)/(2*(K+1)^2)
# Check: the following expression should vanish
  Accuracy+Timeliness+Smoothness+Shift_stopband-MS_error
# Very prototypical: AIC
#  aic<-ifelse(degrees_freedom<K+1&degrees_freedom>1,log(rever)+2*(K-degrees_freedom)/(K)+2*(K-degrees_freedom)*(K-degrees_freedom+1)/(K*(degrees_freedom-1)),NA)

  return(list(Accuracy=Accuracy,Smoothness=Smoothness,Timeliness=Timeliness,MS_error=MS_error))
}















