gtm<-function(y,pr.beta=c(2000,2000,36,12,80,1),pr.intv=rep(NA,6),pr.p0=0.75,pr.dist=NULL,niter=2000,burnin=10,two.sided=TRUE,independent=FALSE,model.strong.ase=TRUE,group.distance=c(1,1,0.5)){

  #Grouped Tissue Model (GTM in the paper)
  #Matti Pirinen 1-Apr-2014, 29-Dec-2014
  #Gibbs sampler to classify tissues into three (or two) groups.
  
  #The default is that we use three groups to specify groups for NOASE, MODASE and SNGASE,
  #G0, NOASE: No allele specific expression (ASE), theta is close to 0.5
  #G1, MODASE: Moderate ASE, theta is different from 0.5 and different from extreme values of 0 and 1
  #G2, SNGASE, Strong ASE, theta is extreme, near 0 or 1
  #Priors for allele frequencies in the groups are given by (truncated) Beta distributions explained below.
  #It is possible to restrict the model to only two groups: NOASE and MODASE.
  
  #The combined configuration is classified into 6 states:
  #C1=NOASE (all tissues in G0),
  #C2=MODASE (all tissues in G1),
  #C3=SNGASE (all tissues in G2)
  #C4=HET0 (heterogeneous with at least one tissue in G0)
  #C5=HET1 (heterogeneous with no tissue in G0)
  #C6=TIS_SPE (tissue specific, one tissue is in different group than all the rest which are in a same group)
  #States C1,...,C5 are non-overlapping and exhaustive, while C6 is a subset of C5.
  #Prior probabilities for states are determined by 'pr.p0' and 'pr.dist' explained below.
  #If only two groups are used, then only combined configurations C1,C2,C4 and C6 are possible.
  
  #INPUT
  #y is m x 2 matrix, of read counts for the two alleles for 'm' tissues (use rownames to identify tissues)
  
  #If model.strong.ase==TRUE then the tissues can come from three groups:
  #theta is the frequency of allele in column 1 of y
  #If two.sided is FALSE then
   #G0: theta~Beta(pr.beta[1],pr.beta[2])*I(pr.intv[1],pr.intv[2]) or if pr.beta[1]==pr.beta[2]==NULL then theta=0.5 (point mass at 0.5)
   #G1: theta~Beta(pr.beta[3],pr.beta[4])*I(pr.intv[3],pr.intv[4])
   #G2: theta~Beta(pr.beta[5],pr.beta[6])*I(pr.intv[5],pr.intv[6])
  #if two.sided is TRUE then
   #G0: theta~0.5*Beta(pr.beta[1],pr.beta[2])*I(pr.intv[1],pr.intv[2])+0.5*Beta(pr.beta[2],pr.beta[1])*I(pr.intv[2],pr.intv[1])
         #or if pr.beta[1]==pr.beta[2]==NULL then theta=0.5
   #G1: theta~0.5*Beta(pr.beta[3],pr.beta[4])*I(pr.intv[3],pr.intv[4])+0.5*Beta(pr.beta[4],pr.beta[3])*I(pr.intv[4],pr.intv[3])
   #G2: theta~0.5*Beta(pr.beta[5],pr.beta[6])*I(pr.intv[5],pr.intv[6])+0.5*Beta(pr.beta[6],pr.beta[5])*I(pr.intv[6],pr.intv[5])
  #If model.strong.ase==FALSE then the tissues come from only two groups: G0 and G1. 
   #In this case 'pr.beta' and 'pr.intv' can be of length 4, or if they are of length 6 then they are truncated to subvector 1:4.
  # 'pr.intv', the three or two intervals (depending on 'model.strong.ase') on which priors for thetas are truncated. If NA, then priors are not truncated.
  # 'independent', if TRUE then each tissue has its own theta, independent of the other tissues (given the priors above)
  #               if FALSE (default) then all tissues in the same group have the same theta
  # 'pr.p0', sum of the prior probability of combined configurations C1,C2, and C3, (see above)
   #if model.strong.ase==TRUE then each of these three configurations has the same prior probability of pr.p0/3 
   #if model.strong.ase==FALSE then configurations C1 and C2 have the same prior probability of pr.p0/2
  # 'pr.dist', vector of prior probability assigned for each set of combined configurations with fixed distance (>0)
   #from homogeneity. These values are Interpreted as relative probabilities and normalised to sum up to (1-pr.p0).
   #Distance of a configuration is the smallest number of tissues whose labels need to be changed 
   #to turn the configuration into a homogeneous one. (E.g. (0,1,1) and (2,1,2) have d=1 and (0,1,2) has d=2.)
   #Each configuration (with distance > 0) will have a prior probability pr.dist[d]/(No. of configurations with distance=d). 
   #Thus two configs with the same distance are equally probable a priori.
   #If model.strong.ase==TRUE then 'pr.dist' should have a length of m-ceiling(m/3).
   #If model.strong.ase==FALSE then 'pr.dist' can have a length of m-ceiling(m/2) or if it has length m-ceiling(m/3) then the last elements are ignored.
   #If pr.dist==NULL, then each distance is given the same prior probability (1-pr.p0)/(m-ceiling(m/3)) or (1-pr.p0)/(m-ceiling(m/2))
   # depending whether 'model.strong.ase' is TRUE or FALSE, respectively.
   #Note that 'pr.p0' is the sum of prior probabilities assigned to the configurations with distance 0 from homogeneity
   #and it is given as a separate parameter, not as part of pr.dist vector.
  # 'niter', the number of Gibbs sampling iterations
  # 'burnin', number of initial iterations discarder, not included in niter (so sampler runs niter+burnin iters in total) 
  #group.distance is a vector of length 3 that tells the distances between groups G0 & G1, G0 & G2 and G1 & G2, respectively
   #these are used to estimate the distances between tissues when model.strong.ase==TRUE.
   #If model.strong.ase==FALSE, then 'group.distance' is ignored and distance between groups G0 and G1 is 1.

  #OUTPUT
  # 'parameters', a list of input parameters
  # 'indiv.posteriors', 3 x m matrix for posterior probabilities of each tissue (col) belonging to each group (rows)
  #row 1 is for G0 (NOASE), row 2 is for G1 (MODASE) and row 3 is for G2 (SNGASE)
  
  # 'distances', an m x m matrix of posterior mean distances between any two tissues
   #when distance between two tissues is as given by 'group.distances'

  # 'top.het.model', the vector of group labels for (a best guess of) the configuration that has maximum likelihood

  # 'log10bfs', log10 of Bayes factors against the null model (=C1, all tissues in NOASE group G0) for the following configurations:
  # 'NOASE', for C1, NOASE, this is always 0
  # 'TOPHET' for configuration 'top.het.model'
  # 'MODASE' for C2, MODASE (all in G1)
  # 'SNGASE' for C3, SNGASE (all in G2)
  # 'HET0' for C4, HET0 (at least one tissue in G0 and one tissue not in G0)
  # 'HET1' for C5, HET1 (none of the tissues in G0 and at least one tissue in G1 and one in G2)
  # 'TIS_SPE' for C6, TIS_SPE (one tissue is different from the others which are in the same group, i.e. configs that have d=1)

  # 'state.posteriors' posterior probabilities for the states C1,...,C6, named as in 'log10bfs'

  
  stopifnot(niter>0 & burnin>0)
  stopifnot(is.logical(model.strong.ase))
  stopifnot(!model.strong.ase | (length(pr.beta)==6 & all(pr.beta[3:6]>0)))
  stopifnot(model.strong.ase | (length(pr.beta) %in% c(4,6) & all(pr.beta[3:4]>0)))
  stopifnot((!is.finite(pr.beta[1]) & !is.finite(pr.beta[2])) | (is.finite(pr.beta[1]) & is.finite(pr.beta[2])))
  if(is.finite(pr.beta[1])) stopifnot(pr.beta[1]>0 & pr.beta[2]>0)
  stopifnot(!model.strong.ase | length(pr.intv)==6)
  stopifnot(model.strong.ase | length(pr.intv)%in% c(4,6))
  stopifnot(is.logical(two.sided) & is.logical(independent))
  stopifnot(ncol(y)==2)
  stopifnot(!model.strong.ase | length(group.distance)==3)
  if(!model.strong.ase) {
    group.distance=rep(1,3)} #d(G0,G1)=1 when only two groups present
  m=nrow(y)
  if(length(rownames(y))==0) rownames(y)=paste("T",1:m,sep="")

  if(m==1) return(gtm.single(y,pr.beta,pr.intv,two.sided,model.strong.ase)) #MP: change this function!
  
  log.prior.dist=logprior.distance(m,p0=pr.p0,p.dist=pr.dist,model.strong.ase=model.strong.ase)
  if(model.strong.ase) log.prior=c(log(pr.p0/3),log.prior.dist$log.prior) #for states with dist==0 prior is pr.p0/3, otherwise from 'logprior.distance'
  if(!model.strong.ase) log.prior=c(log(pr.p0/2),log.prior.dist$log.prior) #for states with dist==0 prior is pr.p0/2, otherwise from 'logprior.distance'
  log.sum.prior.h0=log.prior.dist$log.sum.prior.h0
  log.sum.prior.h1=log.prior.dist$log.sum.prior.h1

  gr=rep(0,m)
  for(i in 1:m){gr[i]=which.max(unlist(gtm.single(matrix(y[i,],nrow=1),pr.beta,pr.intv,two.sided,model.strong.ase)$indiv.posteriors))-1}
  #group indicators, start everything from marginally top state
  #print(gr)
  
  prob=matrix(0,ncol=m,nrow=3) #class probabilities, cols for tissues rows for groups
  loglk=rep(-Inf,3) #log-likelihood + log-prior for three possible states for one tissue in a Gibbs sampling
  dist.01=matrix(0,m,m);dist.02=matrix(0,m,m);dist.12=matrix(0,m,m);
  
  states=rep(NA,(burnin+niter)*m);logprlk=rep(NA,(burnin+niter)*m) #keep track of heterogeneous states and their prior*mlks
                                        #to compute lower bounds for mlks of the heterogeneity states
  in.state=rep(0,5)
  possible.groups=c(0,1,2)
  if(!model.strong.ase) possible.groups=c(0,1)
  for(iter in 1:(burnin+niter)){#start Gibbs

    for(t in 1:m){ #Gibbs update for each tissue separately
      gr.t=gr;#temporary grouping
      for(gr.i in possible.groups){ #tissue t belongs to group gr.i
        gr.t[t]=gr.i;
        loglk[gr.i+1]=logmlk.3.truncated(y=y,gr=gr.t,pr.beta=pr.beta,pr.intv=pr.intv,two.sided=two.sided,independent=independent)+log.prior[1+m-max(as.numeric(table(gr.t)))]
        if(sd(gr.t)>0){
          states[m*(iter-1)+t]=paste(gr.t,collapse="");logprlk[m*(iter-1)+t]=loglk[gr.i+1];}#if heterogeneous state, else leave as NA
      }
      pr=exp(loglk-max(loglk))
      pr=pr/sum(pr)
      gr[t]=sample(c(0,1,2),size=1,prob=pr)#Gibbs sampling from full conditional

      if(iter>burnin) prob[gr[t]+1,t]=prob[gr[t]+1,t]+1
    }
    if(iter>burnin){
      dist.01=dist.01+as.numeric(gr==0)%*%t(as.numeric(gr==1))
      dist.02=dist.02+as.numeric(gr==0)%*%t(as.numeric(gr==2))
      dist.12=dist.12+as.numeric(gr==1)%*%t(as.numeric(gr==2))
      
      if(sd(gr)==0){in.state[gr[1]+1]=in.state[gr[1]+1]+1}
      else{
        if(any(gr==0)) in.state[4]=in.state[4]+1 else in.state[5]=in.state[5]+1}
    }
  } #Gibbs ends
  prob=t(t(prob)/colSums(prob))
  colnames(prob)=rownames(y)
  rownames(prob)=c("NOASE","MODASE","SNGASE")
  dist.01=(dist.01+t(dist.01))/niter
  dist.02=(dist.02+t(dist.02))/niter
  dist.12=(dist.12+t(dist.12))/niter
  distances=group.distance[1]*dist.01+group.distance[2]*dist.02+group.distance[3]*dist.12
  rownames(distances)=rownames(y);colnames(distances)=rownames(y)

  in.state=in.state/niter #currently not returned, as marginal likelihood based posteriors are returned instead
  
  hets=gtm.het(y,prob,pr.beta,pr.intv,two.sided,independent,log.prior)#Go manually through some of the top configurations for heterogeneity states 
  states=c(states,hets$states) #add to the states and logprlk for heterogeneity probability calculations
  logprlk=c(logprlk,hets$logprlk)
  
  top.gr.het=apply(prob,2,which.max)-1 #(best guess for) top model among all heterogeneous states
  if(sd(top.gr.het)==0){#if top model is homogeneous, choose the second best
    prob.2=apply(prob,2,sort)[2,]
    if(sum(prob.2)>0) {
      k=which.max(prob.2)
      top.gr.het[k]=setdiff(order(prob[,k])[2:3]-1,top.gr.het[k])} #this works even when two states have same proba
    else{#whole sampler was in one state, check which tissue should be changed 
      test.states=setdiff(possible.groups,top.gr.het[1]) #checking the other two states
      max.loglk=-Inf;
      for(t in 1:m){ #test each tissue separately
        gr.t=top.gr.het;#temporary grouping
        for(gr.i in test.states){ #tissue t belongs to group gr.i
          gr.t[t]=gr.i;
          loglk=logmlk.3.truncated(y=y,gr=gr.t,pr.beta=pr.beta,pr.intv=pr.intv,two.sided=two.sided,independent=independent) #prior is constant for all these configurations
          if(loglk>max.loglk){max.loglk=loglk;max.t=t;max.state=gr.i}
        }
      }
      top.gr.het[max.t]=max.state
    }
  }
  #print(top.gr)

  if(model.strong.ase){
    top.gr.h1=apply(prob[2:3,],2,which.max) #best guess for top state when none is in state 0
    for(t in 1:m){ #test whether Gibbs is unreliable for each tissue
      if(sum(prob[2:3,t])<1e-2){
        gr.t=top.gr.h1;#temporary grouping
        gr.t[t]=1
        loglk.1=logmlk.3.truncated(y=y,gr=gr.t,pr.beta=pr.beta,pr.intv=pr.intv,two.sided=two.sided,independent=independent)+log.prior[1+m-max(as.numeric(table(gr.t)))]
        gr.t[t]=2
        loglk.2=logmlk.3.truncated(y=y,gr=gr.t,pr.beta=pr.beta,pr.intv=pr.intv,two.sided=two.sided,independent=independent)+log.prior[1+m-max(as.numeric(table(gr.t)))]
        if(loglk.2>loglk.2) top.gr.h1[t]=1 else top.gr.h1[t]=2
      }
    }
    if(sd(top.gr.h1)==0){#if top model is homogeneous, choose the second best
      max.loglk=-Inf;max.t=1;
      for(t in 1:m){
        gr.t=top.gr.h1;gr.t[t]=3-gr.t[t]
        loglk=logmlk.3.truncated(y=y,gr=gr.t,pr.beta=pr.beta,pr.intv=pr.intv,two.sided=two.sided,independent=independent)+log.prior[1+m-max(as.numeric(table(gr.t)))]
        if(loglk>max.loglk){max.loglk=loglk;max.t=t}
      }
      top.gr.h1[max.t]=3-top.gr.h1[t]
    }
  #print(top.gr.h1)
    states=c(states,paste(top.gr.h1,collapse="")) #include in logprlk calculation of het.1 state
    logprlk=c(logprlk,logmlk.3.truncated(y=y,gr=top.gr.h1,pr.beta=pr.beta,pr.intv=pr.intv,two.sided=two.sided,independent=independent)+log.prior[1+m-max(as.numeric(table(top.gr.h1)))])
  }
  
  top.gr.h0=top.gr.het#top state when at least one is in state 0
  if(sum(top.gr.h0==0)==0) {
    max.loglk=-Inf;
    for(t in 1:m){ #test each tissue separately
      gr.t=top.gr.h0;#temporary grouping
      gr.t[t]=0
      loglk=logmlk.3.truncated(y=y,gr=gr.t,pr.beta=pr.beta,pr.intv=pr.intv,two.sided=two.sided,independent=independent)+log.prior[1+m-max(as.numeric(table(gr.t)))]
      if(loglk>max.loglk){max.loglk=loglk;max.t=t;}
    }
    top.gr.h0[max.t]=0
  }
  #print(top.gr.h0)
  states=c(states,paste(top.gr.h0,collapse="")) #include in logprlk calculation of het.0 state
  logprlk=c(logprlk,logmlk.3.truncated(y=y,gr=top.gr.h0,pr.beta=pr.beta,pr.intv=pr.intv,two.sided=two.sided,independent=independent)+log.prior[1+m-max(as.numeric(table(top.gr.h0)))])

  logmlk.state=rep(-Inf,7)
  logmlk.state[1]=logmlk.3.truncated(y=y,gr=rep(0,m),pr.beta=pr.beta,pr.intv=pr.intv,two.sided=two.sided,independent=independent)
  logmlk.state[2]=logmlk.3.truncated(y=y,gr=rep(1,m),pr.beta=pr.beta,pr.intv=pr.intv,two.sided=two.sided,independent=independent)#-loglk.0)/log(10)
  if(model.strong.ase) logmlk.state[3]=logmlk.3.truncated(y=y,gr=rep(2,m),pr.beta=pr.beta,pr.intv=pr.intv,two.sided=two.sided,independent=independent)
  logmlk.state[7]=logmlk.3.truncated(y=y,gr=top.gr.het,pr.beta=pr.beta,pr.intv=pr.intv,two.sided=two.sided,independent=independent)
 
  ind=1;
  if(model.strong.ase){if(m>2) {for.ind=1:m;logmlk.d1=rep(NA,6*m);} else {for.ind=c(1);logmlk.d1=rep(NA,6);}}
  if(!model.strong.ase){if(m>2) {for.ind=1:m;logmlk.d1=rep(NA,2*m);} else {for.ind=c(1);logmlk.d1=rep(NA,2);}}
  for(others in possible.groups){
    for(specific in possible.groups[-(others+1)]){
      for(t in for.ind){
        gr=rep(others,m);gr[t]=specific
        logmlk.d1[ind]=logmlk.3.truncated(y=y,gr=gr,pr.beta=pr.beta,pr.intv=pr.intv,two.sided=two.sided,independent=independent)
        states=c(states,paste(gr,collapse="")) #include in logprlk calculation of het states
        logprlk=c(logprlk,logmlk.d1[ind]+log.prior[2]) #need to add prior for d==1 state here but not in logmlk.d1 above
        #print(c(paste(gr,collapse=""),logmlk.d1[ind]+log.prior[2]))
        ind=ind+1
  }}}
  if(!is.finite(max(logmlk.d1))) {logmlk.state[6]=-Inf} else{
    logmlk.state[6]=max(logmlk.d1)+log(mean(exp(logmlk.d1-max(logmlk.d1))))}#same prior for all d=1 states, use mean as mlk

  #approximate mlk for heterogeneous states from states that Gibbs sampler has visited or that have been separately checked above
  #i.e. assume that all non-visited states have 0 likelihood. Leads to a lower bound.
  keep=!is.na(states)
  states=states[keep];logprlk=logprlk[keep] #remove homogeneous iterations
  keep=!duplicated(states)
  states=states[keep];logprlk=logprlk[keep] #keep only one instance of each state

  #HET0 requires at least one tissue in state 0
  keep=grep("0",states)
  logprlk.k=logprlk[keep];
  #d=rep(NA,length(keep))
  #for(j in 1:length(keep)){
  #  config=as.numeric(unlist(strsplit(states[keep[j]],"")))
  #  d[j]=m-max(as.numeric(table(config)))
  #  print(c(states[keep[j]],logprlk.k[j]))
  #}
  #print(table(d))
  if(!is.finite(max(logprlk.k))) {logmlk.state[4]=-Inf} else{
    logmlk.state[4]=max(logprlk.k)+log(sum(exp(logprlk.k-max(logprlk.k))))-log.sum.prior.h0}

  keep=setdiff(1:length(logprlk),keep)
  logprlk.k=logprlk[keep];
  #d=rep(NA,length(keep))
  #for(j in 1:length(keep)){
  #  config=as.numeric(unlist(strsplit(states[keep[j]],"")))
  #  d[j]=m-max(as.numeric(table(config)))
  #  print(c(states[keep[j]],as.numeric(logprlk.k[j])))
  #}
  #print(table(d))
  #if(!is.finite(max(logprlk.k))) {} #else{
    #logmlk.state[5]=max(logprlk.k)+log(sum(exp(logprlk.k-max(logprlk.k))))-log.sum.prior.h1}
  logmlk.state[5]=-Inf
  prob=data.frame(prob,row.names=c("NOASE","MODASE","SNGASE"));names(prob)=rownames(y)
  logmlk.0=logmlk.state[min(which(is.finite(logmlk.state)))] #reference state, preferably NOASE if it is not -Inf
  log10bfs=data.frame(matrix((logmlk.state-logmlk.0)/log(10),nrow=1))
  names(log10bfs)=c("NOASE","MODASE","SNGASE","HET0","HET1","TIS_SPE","TOPHET")
  if(model.strong.ase) log10.prior.states=log10(c(rep(pr.p0/3,3),exp(c(log.sum.prior.h0,log.sum.prior.h1))))
  if(!model.strong.ase) log10.prior.states=log10(c(rep(pr.p0/2,2),0,exp(c(log.sum.prior.h0,log.sum.prior.h1))))
  if(model.strong.ase) {if(m>2) log10.prior.d1=(log.prior[2]+log(6*m))/log(10) else log10.prior.d1=(log.prior[2]+log(6))/log(10)}
  if(!model.strong.ase) {if(m>2) log10.prior.d1=(log.prior[2]+log(2*m))/log(10) else log10.prior.d1=(log.prior[2]+log(2))/log(10)}
 
  posteriors=log10.prior.states+log10bfs[1,1:5]
  posteriors=10^(posteriors-max(posteriors))
  posteriors=data.frame(matrix(posteriors/sum(posteriors),nrow=1))
  i=which.max(unlist(posteriors[1,]))
  posteriors=cbind(posteriors,min(c(1,10^(log10(unlist(posteriors[1,i]))+log10bfs[1,6]-log10bfs[1,i]+log10.prior.d1-log10.prior.states[i])))) #p(d=1|data)=p(i|data)* p(d=1)/p(i) * p(data|d=1)/p(data|i)
  names(posteriors)=c("NOASE","MODASE","SNGASE","HET0","HET1","TIS_SPE")

  parameters=list(pr.beta=pr.beta,pr.intv=pr.intv,pr.p0=pr.p0,pr.dist=pr.dist,niter=niter,two.sided=two.sided,independent=independent,model.strong.ase=model.strong.ase,group.distance=group.distance)
  
  return(list(parameters=parameters,indiv.posteriors=prob,distances=distances,top.het.model=top.gr.het,log10bfs=log10bfs,state.posteriors=posteriors))
  
}


gtm.star<-function(y.list,pr.beta=c(2000,2000,36,12,80,1),pr.intv=rep(NA,6),pr.pi=rep(1,5),niter=2000,burnin=10,two.sided=TRUE,independent=FALSE,model.strong.ase=TRUE,group.distance=c(1,1,0.5)){

  #Multi-locus Grouped Tissue Model (GTM* in the paper)
  #Matti Pirinen 01-Apr-2014, 29-Dec-2014
  #Hierarchical model where each variant belongs to one of five states and each tissue within a variant belongs to one of three (or two) groups
  
  #We use three groups to specify groups for NOASE, MODASE and SNGASE,
  #G0, NOASE: No allele specific expression (ASE), theta is close to 0.5
  #G1, MODASE: Moderate ASE, theta is different from 0.5 and different from extreme values of 0 and 1
  #G2, SNGASE, Strong ASE, theta is extreme, near 0 or 1
  #Priors for allele frequencies in the groups are given by Beta distributions explained below.
  #It is possible to restrict the model to only two groups: NOASE and MODASE.

  #The combined configuration is classified into 6 states:
  #C1=NOASE (all tissues in G0),
  #C2=MODASE (all tissues in G1),
  #C3=SNGASE (all tissues in G2)
  #C4=HET0 (heterogeneous with at least one tissue in G0)
  #C5=HET1 (heterogeneous with no tissue in G0)
  #C6=TIS_SPE (tissue specific, one tissue is in different group than all the rest which are in a same group)
  #States C1,...,C5 are non-overlapping and exhaustive, while C6 is a subset of C5.
  #Prior probabilities for states C1,...,C5 are ~ Dirichlet(pr.pi)
  #If only two groups are used, then only combined configurations C1,C2,C4 and C6 are possible.

  #INPUT
  #y.list is a list of matrices where
  #y.list[[s]] is an m[s] x 2 matrix, of read counts for the two alleles for 'm[s]' tissues
  #these matrices MUST have rownames so that the tissues can be matched across the variants
  
  #If model.strong.ase==TRUE then the tissues can come from three groups:
  #theta is the frequency of allele in column 1 of y
  #if two.sided is FALSE then
   #G0: theta~Beta(pr.beta[1],pr.beta[2])*I(pr.intv[1],pr.intv[2]) or if pr.beta[1]==pr.beta[2]==NULL then theta=0.5 (point mass at 0.5)
   #G1: theta~Beta(pr.beta[3],pr.beta[4])*I(pr.intv[3],pr.intv[4])
   #G2: theta~Beta(pr.beta[5],pr.beta[6])*I(pr.intv[5],pr.intv[6])
  #if two.sided is TRUE then
   #G0: theta~0.5*Beta(pr.beta[1],pr.beta[2])*I(pr.intv[1],pr.intv[2])+0.5*Beta(pr.beta[2],pr.beta[1])*I(pr.intv[2],pr.intv[1])
         #or if pr.beta[1]==pr.beta[2]==NULL then theta=0.5
   #G1: theta~0.5*Beta(pr.beta[3],pr.beta[4])*I(pr.intv[3],pr.intv[4])+0.5*Beta(pr.beta[4],pr.beta[3])*I(pr.intv[4],pr.intv[3])
   #G2: theta~0.5*Beta(pr.beta[5],pr.beta[6])*I(pr.intv[5],pr.intv[6])+0.5*Beta(pr.beta[6],pr.beta[5])*I(pr.intv[6],pr.intv[5])
  #If model.strong.ase==FALSE then the tissues come from only two groups: G0 and G1. 
   #In this case 'pr.beta' and 'pr.intv' can be of length 4, or if they are of length 6 then they are truncated to subvector 1:4.
  # 'pr.intv', the three intervals on which priors for thetas are truncated. If NA, then priors are not truncated.
  # 'independent', if TRUE then each tissue has its own theta, independent of the other tissues (given the priors above)
  #               if FALSE (default) then all tissues in the same group have the same theta
  # 'niter', the number of MCMC iterations
  # 'burnin', number of initial iterations discarder, not included in niter (so sampler runs niter+burnin iters in total) 
  # 'pr.pi', parameters for Dirichlet prior on probability of each of the five states
   #OR if length(pr.pi)=4 then HET0 and HET1 states are combined to a single HET state
   #If model.strong.ase==TRUE then only values 1,2 and 4 of pr.pi are used. In this case it is also possible to give just three values fro pr.pi.
  #prior probability of each heterogeneous configuration depends on pi[4] or pi[5] and the distance of the configuration.
   #Distance of a configuration is the smallest number of tissues whose labels need to be changed 
   #to turn the configuration into a homogeneous one. (E.g. (0,1,1) and (2,1,2) have d=1 and (0,1,2) has d=2.)
   #If model.strong.ase==TRUE then maximum.distance=m-ceiling(m/3), else maximum.distance=floor(m/2).
   #if model.strong.ase=FALSE or length(pr.pi)==4 then
    #Each configuration (with distance d > 0) will have a prior probability pi[4]/maximum.distance/(No. of configurations with distance=d).
    #Thus two configs with the same distance are equally probable a priori.
   #if length(pr.pi)==5 then
    #there are separate priors for HET0 and HET1 configurations, pi[4](m-ceiling(m/3))/(No. of tissues with distance=d)
    #and pi[5]/(floor(m/2))/(No. of tissues with distance=d), respectively
  #group.distance is a vector of length 3 that tells the distances between groups 0-1, 0-2 and 1-2, respectively
  #these are used to estimate the distances between tissues
  
  
  #OUTPUT
  #list 'indiv.posteriors', variant specific 3 x m matrices for posterior probabilities of each tissue (col) belonging to each group (rows)
  #row 1 is for G0 (NOASE), row 2 is for G1 (MODASE) and row 3 is for G2 (SNGASE)
  
  #list 'state.posteriors' variant specific posterior probabilities for the states C1,...,C6

  # 'prop.posteriors' 5 x 3 matrix. col1 posterior for pi[1,...,5], col2=lower 95% point, col3=upper 95% point

  #distances, a matrix of nt x nt where nt is the total number of tissues across all variants s=1,...,length(y.list)
   #each value is an estimated distance between two tissues, where 0 means that tissues are always in the same group
   #if there are no common variants for two tissues, then distance is NA

  #distances.se, a matrix of nt x nt where nt is the total number of tissues across all variants s=1,...,length(y.list)
   #each value is a standard error for distance estimate according to multinomial sampling model where each variant
   #corresponds to one observation
  
  stopifnot(niter>0 & burnin>0)
  stopifnot(is.logical(model.strong.ase))
  stopifnot(!model.strong.ase | (length(pr.beta)==6 & all(pr.beta[3:6]>0)))
  stopifnot(model.strong.ase | (length(pr.beta) %in% c(4,6) & all(pr.beta[3:4]>0)))
  stopifnot((!is.finite(pr.beta[1]) & !is.finite(pr.beta[2])) | (is.finite(pr.beta[1]) & is.finite(pr.beta[2])))
  if(is.finite(pr.beta[1])) stopifnot(pr.beta[1]>0 & pr.beta[2]>0)
  stopifnot(!model.strong.ase | length(pr.intv)==6)
  stopifnot(model.strong.ase | length(pr.intv)%in% c(4,6))
 
  stopifnot(is.logical(two.sided) & is.logical(independent))
  stopifnot(all(unlist(lapply(y.list,ncol))==2))
  if(!model.strong.ase & length(pr.pi)==3) pr.pi[4]=pr.pi[3] #copy from 3 to 4 and use only 1,2 and 4 below
  stopifnot(all(pr.pi>0) & length(pr.pi) %in% c(4,5))
  stopifnot(!model.strong.ase | length(group.distance)==3)
  if(!model.strong.ase) {
    group.distance=rep(1,3)} #d(G0,G1)=1 when only two groups present

  nsets=length(y.list)
  m=unlist(lapply(y.list,nrow))
  stopifnot(all(m>1))

  log.prior=list();log.prior.h0=list();log.prior.h1=list();gr=list();prob=list();state.prob=list();tissues=c();
  dist.0=list();dist.01=list();dist.02=list();dist.12=list();
  for(i in 1:nsets){
    stopifnot(rownames(y.list[[i]])!=NULL) #tissues must have names
    log.prior.dist=logprior.distance(m[i],p0=0.75,p.dist=NULL,model.strong.ase=model.strong.ase)
    log.prior[[i]]=log(4)+log.prior.dist$log.prior #by multiplying this prior by p.het, we get prior for each het distance, 4 is because (1-0.75)=1/4
    log.prior.h0[[i]]=-log.prior.dist$log.sum.prior.h0+log.prior.dist$log.prior.h0 #by multiplying this prior by p.h0, we get prior for each het0 distance
    log.prior.h1[[i]]=-log.prior.dist$log.sum.prior.h1+log.prior.dist$log.prior.h1 #by multiplying this prior by p.h1, we get prior for each het1 distance
    gr[[i]]=rep(0,m[i]) #group indicators, start everything from group 0
    for(s in 1:m[i]){
        gr[[i]][s]=which.max(unlist(gtm.single(matrix(y.list[[i]][s,],nrow=1),pr.beta,pr.intv,two.sided,model.strong.ase)$indiv.posteriors))-1
    }
    prob[[i]]=matrix(0,ncol=m[i],nrow=3) #class probabilities, cols for tissues rows for groups
    state.prob[[i]]=rep(0,6)
    dist.0[[i]]=matrix(0,m[i],m[i]);dist.01[[i]]=matrix(0,m[i],m[i]);dist.02[[i]]=matrix(0,m[i],m[i]);dist.12[[i]]=matrix(0,m[i],m[i])
    tissues=c(tissues,rownames(y.list[[i]]))
  }
  tissues=sort(unique(tissues))
  loglk=rep(-Inf,3) #log-likelihood + log-prior for three possible states for one tissue in a Gibbs sampling
  possible.groups=c(0,1,2)
  if(!model.strong.ase){
    possible.groups=c(0,1)
    pr.pi[3]=0
    pr.pi=pr.pi[1:4]
  }
  pi=pr.pi/sum(pr.pi)
  npi=length(pi)
  logpi=log(pi)
  pi.res=matrix(NA,ncol=npi,nrow=niter)
  
  for(iter in 1:(burnin+niter)){#start Gibbs
    print(str_c("iter = ", iter))
    nstate=rep(0,npi) #how many in each state 1=NOASE,2=MODASE,3=SNGASE,4=HET (or 4=HET0,5=HET1)
    for(s in 1:nsets){
      for(t in 1:m[s]){ #Gibbs update for each tissue separately
        gr.t=gr[[s]];#temporary grouping
        for(gr.i in possible.groups){ #tissue t belongs to group gr.i
          gr.t[t]=gr.i;
          loglk[gr.i+1]=logmlk.3.truncated(y=y.list[[s]],gr=gr.t,pr.beta=pr.beta,pr.intv=pr.intv,two.sided=two.sided,independent=independent)
          dist=m[s]-max(as.numeric(table(gr.t)))
          if(dist>0) {
            if(npi==4) {loglk[gr.i+1]=loglk[gr.i+1]+logpi[4]+log.prior[[s]][dist]} else{
              if(any(gr.t==0)) {loglk[gr.i+1]=loglk[gr.i+1]+logpi[4]+log.prior.h0[[s]][dist]}else{
                loglk[gr.i+1]=loglk[gr.i+1]+logpi[5]+log.prior.h1[[s]][dist]}
            }
          }
          else{loglk[gr.i+1]=loglk[gr.i+1]+logpi[1+gr.i]}          
        }
        pr=exp(loglk-max(loglk))
        pr=pr/sum(pr)
        gr[[s]][t]=sample(c(0,1,2),size=1,prob=pr)#Gibbs sampling from full conditional
        if(iter>burnin) prob[[s]][gr[[s]][t]+1,t]=prob[[s]][gr[[s]][t]+1,t]+1
      }
      if(sd(gr[[s]])==0){ state=1+gr[[s]][1]}else{
        if(any(gr[[s]]==0)) state=4 else state=5}
      if(npi==4) nstate[c(1:4,4)[state]]=1+nstate[c(1:4,4)[state]] else nstate[state]=nstate[state]+1
      if(iter>burnin){
        state.prob[[s]][state]=1+state.prob[[s]][state]
        state.prob[[s]][6]=state.prob[[s]][6]+((m[s]-max(as.numeric(table(gr[[s]]))))==1) #TIS_SPE
        dist.01[[s]]=dist.01[[s]]+as.numeric(gr[[s]]==0)%*%t(as.numeric(gr[[s]]==1))+as.numeric(gr[[s]]==1)%*%t(as.numeric(gr[[s]]==0))
        dist.02[[s]]=dist.02[[s]]+as.numeric(gr[[s]]==0)%*%t(as.numeric(gr[[s]]==2))+as.numeric(gr[[s]]==2)%*%t(as.numeric(gr[[s]]==0))
        dist.12[[s]]=dist.12[[s]]+as.numeric(gr[[s]]==1)%*%t(as.numeric(gr[[s]]==2))+as.numeric(gr[[s]]==2)%*%t(as.numeric(gr[[s]]==1))
      }
    }
    #update pi, start only after burnin
    if(iter>burnin) {
      if(model.strong.ase) for(i in 1:npi) pi[i]=rgamma(1,shape=pr.pi[i]+nstate[i],scale=1)
      if(!model.strong.ase) for(i in c(1,2,4)) pi[i]=rgamma(1,shape=pr.pi[i]+nstate[i],scale=1)
    }
    pi=pi/sum(pi)
    logpi=log(pi)
    if(iter>burnin) pi.res[iter-burnin,]=pi
  } #Gibbs ends
  
  nallt=length(tissues)
  mean.d=matrix(0,nallt,nallt);var.d=matrix(0,nallt,nallt);nsamples=matrix(0,nallt,nallt);
  for(s in 1:nsets){
    prob[[s]]=t(t(prob[[s]])/colSums(prob[[s]]))
    rownames(prob[[s]])=c("NOASE","MODASE","SNGASE")
    colnames(prob[[s]])=rownames(y.list[[s]])
    state.prob[[s]]=state.prob[[s]]/niter
    names(state.prob[[s]])=c("NOASE","MODASE","SNGASE","HET0","HET1","TIS_SPE")
    dist.01[[s]]=dist.01[[s]]/niter
    dist.02[[s]]=dist.02[[s]]/niter
    dist.12[[s]]=dist.12[[s]]/niter
    dist.0[[s]]=1-dist.01[[s]]-dist.02[[s]]-dist.12[[s]]
    nams=rownames(y.list[[s]])
    for(i in 2:m[s]){
      i.ind=which(tissues == nams[i]) 
      for(j in 1:(i-1)){
        j.ind=which(tissues == nams[j])
        p.vec=c(dist.0[[s]][i,j],dist.01[[s]][i,j],dist.02[[s]][i,j],dist.12[[s]][i,j])
        mean.d[i.ind,j.ind]=mean.d[i.ind,j.ind]+t(p.vec)%*%c(0,group.distance)
        v.mat=-p.vec%*%t(p.vec)
        diag(v.mat)=p.vec*(1-p.vec)
        var.d[i.ind,j.ind]=var.d[i.ind,j.ind]+t(c(0,group.distance)) %*% v.mat %*% c(0,group.distance)
        nsamples[i.ind,j.ind]=nsamples[i.ind,j.ind]+1
      }
    }
  }
  nsamples=nsamples+t(nsamples)
  mean.d=(mean.d+t(mean.d))/nsamples;var.d=(var.d+t(var.d))/nsamples;
  mean.d[!is.finite(mean.d)]=NA;var.d[!is.finite(mean.d)]=NA;
  diag(mean.d)=0;diag(var.d)=0;
  rownames(mean.d)=tissues;colnames(mean.d)=tissues;
  rownames(var.d)=tissues;colnames(var.d)=tissues;
  
  hm.res=matrix(NA,nrow=npi,ncol=3)
  hm.res[,1]=as.numeric(colSums(pi.res)/niter)
  hm.res[,2]=as.numeric(apply(pi.res,2,function(x) quantile(x,0.025)))
  hm.res[,3]=as.numeric(apply(pi.res,2,function(x) quantile(x,0.975)))
  if(npi==5) rownames(hm.res)=c("NOASE","MODASE","SNGASE","HET0","HET1")
  if(npi==4) rownames(hm.res)=c("NOASE","MODASE","SNGASE","HET")
  colnames(hm.res)=c("est","low95","up95")

  parameters=list(pr.beta=pr.beta,pr.intv=pr.intv,niter=niter,two.sided=two.sided,independent=independent,model.strong.ase=model.strong.ase,group.distance=group.distance)
  
  return(list(parameters=parameters,indiv.posteriors=prob,state.posteriors=state.prob,distances=mean.d,distances.se=sqrt(var.d),prop.posteriors=hm.res))
}

gtm.single<-function(y,pr.beta,pr.intv,two.sided,model.strong.ase){
  #returns group probabilities for a single tissue
  stopifnot(ncol(y)==2 & nrow(y)==1)
  logmlk=rep(-Inf,3)
  #prior for each state is the same
  possible.groups=0:2
  if(!model.strong.ase) possible.groups=0:1
  for(gr in possible.groups){
    logmlk[gr+1]=logmlk.3.truncated(y=y,gr=gr,pr.beta=pr.beta,pr.intv=pr.intv,two.sided=two.sided,independent=TRUE)}
  i=min(which(is.finite(logmlk)))
  log10bfs=data.frame(matrix((logmlk-logmlk[i])/log(10),nrow=1))
  names(log10bfs)=c("NOASE","MODASE","SNGASE")
  posteriors=log10bfs
  posteriors=10^(log10bfs-max(log10bfs))
  posteriors=data.frame(matrix(posteriors/sum(posteriors),nrow=1))
  names(posteriors)=c("NOASE","MODASE","SNGASE")
  prob=matrix(posteriors[1,],ncol=1);rownames(prob)=c("NOASE","MODASE","SNGASE");colnames(prob)=rownames(y)
  parameters=list(pr.beta=pr.beta,pr.intv=pr.intv,two.sided=two.sided,model.strong.ase=model.strong.ase)
  return(list(parameters=parameters,indiv.posteriors=prob,log10bfs=log10bfs,state.posteriors=posteriors))
  
}

gtm.het<-function(y,prob,pr.beta,pr.intv,two.sided,independent,log.prior){
  #To check some of the most probable heterogeneity configurations based on individual probabilities
  #Taken these into account in marg likelihood calculation will improve accuracy of het probabilities
  
  n.max.tis=5 #max n of tissues that can be in two groups, at most 2^n.max.tis different configs will be checked
  uncertainty.threshold=0.2 #if the second best group has proba over this, then both groups are considered  
  nt=ncol(prob)
  stopifnot(nrow(prob)==3)
  stopifnot(length(log.prior) %in% c(1+nt-ceiling(nt/3),1+floor(nt/2)))
  
  prob.2=apply(prob,2,sort)[2,]
  het.tis=rev(order(prob.2))
  k=max(c(2,min(c(n.max.tis,sum(prob.2[het.tis]>uncertainty.threshold)))))
  het.tis=het.tis[1:k]
  #go through all 2^k possible heterogeneity states for het.tis
    
  two.groups=apply(prob[,het.tis],2,order)[2:3,]-1
  top.groups=apply(prob,2,order)[3,]-1
  states=rep(NA,2^k)
  logprlk=rep(NA,2^k)
  
  j=1
  i1<-2^(0:(k-1))
  i2<-2*i1
  for(top.set in 0:(2^k-1)){#any subset can be in their top groups
    gr=as.numeric(top.groups)
    ind=(top.set %% i2)>=i1
    gr[het.tis[ind]]=two.groups[2,ind] #these are in their top state
    gr[het.tis[!ind]]=two.groups[1,!ind] #these are in their 2nd best state
    #print(gr)
    if(sd(gr)>0){ #leave homogeneous configurations to NA
      states[j]=paste(gr,collapse="")
      logprlk[j]=logmlk.3.truncated(y,gr,pr.beta,pr.intv,two.sided,independent)+log.prior[1+nt-max(as.numeric(table(gr)))]
    }
    j=j+1
  } 
  return(list(states=states,logprlk=logprlk))
}
   
      
logmlk.3.truncated<-function(y,gr,pr.beta,pr.intv,two.sided,independent){
  #returns marginal loglikelihood,
  #for input parameters see 'gtm'
  #allows for truncated priors and independence across tissues

  #stopifnot(length(pr.intv)==6)
  #stopifnot(length(pr.beta)==6)
  stopifnot(all(gr %in% c(0,1,2)))
  
  loglk=0
  if(!is.finite(pr.beta[1]) | !is.finite(pr.beta[2])){
    for.inds=c(1,2)
    k=which(gr==0)
    if(length(k)>0){loglk=loglk+sum(y[k,])*log(0.5)}
  }else{for.inds=c(0,1,2) }
  
  for(gr.i in for.inds){ #either gr.i in c(1,2) or gr.i in c(0,1,2)
    k=which(gr==gr.i)
    if(length(k)>0){
      aa=pr.beta[gr.i*2+1]
      bb=pr.beta[gr.i*2+2]

      lbeta.pr=lbeta(aa,bb); #this is the same as lbeta(bb,aa)
      x0=pr.intv[2*gr.i+1];x1=pr.intv[2*gr.i+2]
      if(!is.finite(x0) | x0<0) x0=0; if(!is.finite(x1) | x1>1) x1=1; 
      stopifnot(x0<x1)
      if(independent) yy=matrix(y[k,],ncol=2,byrow=FALSE) else yy=matrix(c(sum(y[k,1]),sum(y[k,2])),nrow=1)

      for(i in 1:nrow(yy)){
      #one sided prior
        lp1=pbeta(x1,aa,bb,log=TRUE)
        lp0=pbeta(x0,aa,bb,log=TRUE)
        lp.pr.one=lp1+log(1-exp(lp0-lp1)) #=log(p1-p0)
      #one sided posterior:
        lp1=pbeta(x1,aa+yy[i,1],bb+yy[i,2],log=TRUE)
        lp0=pbeta(x0,aa+yy[i,1],bb+yy[i,2],log=TRUE)
        lp.one=lp1+log(1-exp(lp0-lp1)) #=log(p1-p0)
   
        if(!two.sided){
          loglk=loglk+lbeta(aa+yy[i,1],bb+yy[i,2])+lp.one-lbeta.pr-lp.pr.one}
        if(two.sided){
          u1=lbeta(aa+yy[i,1],bb+yy[i,2])+lp.one-lbeta.pr-lp.pr.one
        #second side of the two-sided prior, interval is (1-x1,1-x0)
          lp1=pbeta(1-x0,bb,aa,log=TRUE)
          lp0=pbeta(1-x1,bb,aa,log=TRUE)
          lp.pr.two=lp1+log(1-exp(lp0-lp1)) #=log(p1-p0)
        #posterior:
          lp1=pbeta(1-x0,bb+yy[i,1],aa+yy[i,2],log=TRUE)
          lp0=pbeta(1-x1,bb+yy[i,1],aa+yy[i,2],log=TRUE)
          lp.two=lp1+log(1-exp(lp0-lp1)) #=log(p1-p0)
          u2=lbeta(bb+yy[i,1],aa+yy[i,2])+lp.two-lbeta.pr-lp.pr.two 
          u=max(c(u1,u2))
          if(!is.finite(u)) loglk=-Inf else loglk=loglk+u+log(0.5*(exp(u1-u)+exp(u2-u)))
        }
      }
    }
  }
  return(loglk)
}


count.states<-function(m,model.strong.ase=TRUE){
  #MP 24-Feb-2014, 29-Dec-2014
  #Counts how many heterogeneous states there are
  #as a function of distance from homogeneous states.
  #Distance is the smallest number of changes that turns the state into one of the homogeneous states.
  #Thus maximum distance is m-ceiling(m/3) when model.strong.ase==TRUE, and floor(m/2) when model.strong.ase==FALSE
  # minimum is 1 (for a heterogeneous state).

  #INPUT
  # 'm' the number of tissues.
  #model.strong.ase, if TRUE then there are three groups, if FALSE then only two groups
  
  #OUTPUT
  # 'logn.all' log of the number of all heterogeneous states. 
  # 'logn.0' log of the number of heterogeneous states belonging to category HET0, i.e. at least one 0.
  # 'logn.1' log of the number of heterogeneous states belonging to category HET1, i.e. none in 0.

  if(!model.strong.ase){ #only two groups, all heterogeneous belong to HET0 and non to HET1
    logn.all=rep(NA,floor(m/2))
    if(m>2){for(k in 1:floor((m-1)/2)) {logn.all[k]=log(2)+lchoose(m,k)}}
    if(m%%2==0) logn.all[m/2]=lchoose(m,m/2)
    logn.0=logn.all
    logn.1=log(rep(0,floor(m/2)))
  }
    
  if(model.strong.ase){
    dvals.all=seq(1,m-ceiling(m/3))
    logn.all=rep(0,length(dvals.all))
    for(k in dvals.all){
      i.list=max(c(0,2*k-m)):min(c(k,m-k))
      i.list=i.list[(k-i.list>=i.list)]
      logn.i=rep(NA,length(i.list))
      j=1
      for(i in i.list){
        logn.i[j]=lfactorial(m)-lfactorial(m-k)-lfactorial(i)-lfactorial(k-i)
        vs=c(m-k,i,k-i)
        logn.i[j]=logn.i[j]+log(c(1,3,6)[length(unique(vs))])
        j=j+1
      }
      logn.all[k]=log(sum(exp(logn.i)))
    }
    #print(paste("Count:",sum(exp(logn.all))," exact=",(3^m-3)))

    dvals.1=seq(1,floor(m/2))
    logn.1=rep(0,length(dvals.1))
    for(k in dvals.1){
      logn.1[k]=lfactorial(m)-lfactorial(m-k)-lfactorial(k)
      if(k<(m/2)) logn.1[k]=logn.1[k]+log(2)
    }
    #print(paste("Count:",sum(exp(logn.1))," exact=",(2^m-2)))
  
    logn.0=log(exp(logn.all)-c(exp(logn.1),rep(0,length(logn.all)-length(logn.1))))
  }

  return(list(logn.0=logn.0,logn.1=logn.1,logn.all=logn.all))
}


logprior.distance<-function(m,p0=0.75,p.dist=NULL,model.strong.ase=TRUE){
  #Counts the prior probability for each heterogeneous state as a function of its distance. 
  #Distance is the smallest number of changes that turns the state into one of the homogeneous states
  #if model.strong.ase==TRUE then the maximum distance is m-ceiling(m/3),
  #if model.strong.ase==FALSE then the maximum distance is floor(m/2),
  #minimum distance is 1 for a heterogeneous configuration
  #(the three homogeneous configurations have distance 0)
  #INPUT
  # 'm' the number of tissues
  # 'p0' the joint prior probability of the 3 homogeneous states
  # 'p.dist' either a vector of length 'max.dist' of total probabilities of each
  #         set of states for distances 1,...,'max.dist'
  #         Where max.dist==m-ceiling(m/3) if model.strong.ase==TRUE and max.dist==floor(m/2) if model.strong.ase==FALSE.
  #         Interpreted as relative to each other so that after renormalisation and scaling sum to (1-p0).
  #OR, if p.dist == NULL,
  #         then p.dist will be set uniform =(1-p0)/max.dist over the distance.
  #OUTPUT
  # 'log.prior' is the log of prior probability of each STATE as a function of distance 1,...,max.dist
  #            results from dividing p.dist by the number of states in each distance category
  # 'log.sum.prior.h0' is the log of the sum of priors over all heterog states that have at least one 0
  # 'log.sum.prior.h1' is the log of the sum of priors over all heterog states that have no 0
  
  stopifnot( p0<=1 & p0>=0 )
  counts=count.states(m,model.strong.ase)
  if(model.strong.ase){ndist=m-ceiling(m/3)}
  if(!model.strong.ase){ndist=floor(m/2)}
  if(is.null(p.dist)){
    p.dist=rep((1-p0)/ndist,ndist)
  }else{
    stopifnot(length(p.dist)==ndist)
    stopifnot(all(p.dist >= 0))
    p.dist=(1-p0)*p.dist/sum(p.dist) #renormalise
  }
  log.prior=log(p.dist)-counts$logn.all #log prior for a het config as a function of distance when mass distributed among all het configs.
  log.sum.prior.h0=log(sum(exp(counts$logn.0+log.prior)))
  log.prior.h0=log.sum.prior.h0-log(ndist)-counts$logn.0 #log prior for a HET0 config as function of distance when distributed only among HET0
  log.sum.prior.h1=log(sum(exp(counts$logn.1+log.prior[1:length(counts$logn.1)])))
  log.prior.h1=log.sum.prior.h1-log(length(counts$logn.1))-counts$logn.1 #log prior for a HET1 config as function of distance when distributed only among HET1

  if(!is.finite(log.sum.prior.h1)) log.prior.h1=rep(-Inf,length(log.prior.h1))
  
  return(list(log.prior=log.prior,log.sum.prior.h0=log.sum.prior.h0,log.sum.prior.h1=log.sum.prior.h1,log.prior.h0=log.prior.h0,log.prior.h1=log.prior.h1))
}


plot.gtm<-function(res,y,title.text=NULL,print.counts=TRUE,cex.c=1.0){
  #INPUT
  # 'res' list from function 'gtm'
  # 'y' mx2 matrix of allele read counts, col 1 for reference allele, col 2 for non ref allele allele, uses rownames as labels if present
  # 'title.text' if present uses as a title in the plot
  #print.counts if TRUE, prints out the read counts in the form NONREF/ALL, if FALSE does not print counts 
  
  #OUTPUT
  #a plot
  
  m=length(y[,1])
  if(length(rownames(y))>0) tissues=rownames(y) else tissues=paste("T",1:m,sep="")
  
  layout(matrix(c(3,2,1),ncol=1,nrow=3,byrow=TRUE),height=c(1.3,2,1),width=c(4))
  
  par(mar=c(6,5,1,1))
  b.plot=barplot(as.matrix(res$indiv.posteriors),beside=FALSE,col=c("white","grey","black"),ylab="INDIVIDUAL PROB",cex.lab=1.2*cex.c,cex.axis=1.3*cex.c,cex.main=1.3*cex.c,xaxt="n",width=1,yaxt="n",xlim=c(0.2+(m-2)*0.0375,m*1.16))
  text(b.plot,-0.35,tissues,srt=45,xpd=TRUE,cex=1.1*cex.c)
  axis(2,at=c(0,0.25,0.5,0.75,1),labels=c(0,"",0.5,"",1),cex.axis=1.3*cex.c)

  par(mar=c(1,5,1,1))
  plot(0,0,col="white",xlim=c(0,1+max(b.plot)),ylim=c(0,1),xlab="",xaxt="n",ylab="NON REF FREQ",main="",xaxs="i",yaxs="i",cex.lab=1.2*cex.c,cex.axis=1.3*cex.c,cex.main=1.3*cex.c)    
  abline(h=0.5,lty=2)
  coeff=1.05
  if(m>8) coeff=1.02
  if(m>20) coeff=1.01
  if(m>40) coeff=1
  for(t in 1:m){
    x.coord=b.plot[t]*coeff
    post.a=0.5+y[t,1];post.b=0.5+y[t,2]
    freq=y[t,2]/sum(y[t,])
    arrows(x.coord,qbeta(0.025,post.b,post.a),x.coord,qbeta(0.975,post.b,post.a),code=0,angle=90,lwd=1.5*cex.c)
    points(x.coord,freq,pch=19,cex=1.5*cex.c)
    if(print.counts) mtext(paste(y[t,2],"/",sum(y[t,]),sep=""),1,line=0.6,at=x.coord,cex=cex.c)
  }

  if(length(title.text) > 0) par(mar=c(5,7,4,2)) else  par(mar=c(5,7,1,2))
  cols=c("white","grey","black","cyan","violet","springgreen3")
  if(m==1) cols=cols[1:3]
  b.plot=barplot(rev(as.matrix(res$state.posteriors)),horiz=TRUE,beside=TRUE,col=rev(cols),ylab="",cex.lab=1.2*cex.c,cex.axis=1.3*cex.c,cex.main=2*cex.c,xlab="STATE PROBABILITIES",yaxt="n",xlim=c(0,1),main=title.text)  
  text(-0.1,b.plot,rev(names(res$state.posteriors)),srt=0,xpd=TRUE,cex=1.2*cex.c)
  
}



