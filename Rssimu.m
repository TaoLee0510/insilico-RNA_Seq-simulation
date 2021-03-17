function Rssimu(K_data, r_data, n, sampleing_number, mix_ratio, cylce_number)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Rssimu: RNA_Seq simulation%%%%%%%%%%%%%%%%%%%%%%%%%
%K_data : k data
%r_data : r data
%n : numbers of genes (20531)
%sampleing_number
%mix_ratio : ratio of K over r
%cylce_number
%%%%%%%%%%%%%%%%%%%%%
file_K=fopen(K_data);
file_r=fopen(r_data);
n=str2double(n);
m=str2double(sampleing_number);
MixNumber=str2double(mix_ratio);
b=str2double(cylce_number);
xx=10*MixNumber;
yy=10*(1-MixNumber);
dirname=['Kr_' num2str(xx) '_' num2str(yy)];
D=['mkdir ' dirname];
system(D);
CK_all=cell(n+1,3);
Cr_all=cell(n+1,3);
CK=zeros(n,3);
Cr=zeros(n,3);
    for a1=1:1:n+1
        for b1=1:1:3
            CK_all{a1,b1}=fscanf(file_K,'%s\t',[1,1]);
        end
    end
    for c1=1:1:n+1
        for d1=1:1:3
            Cr_all{c1,d1}=fscanf(file_r,'%s\t',[1,1]);
        end
    end   
    
for e1=2:1:n
    for f1=2:1:3
        CK(e1-1,f1)=str2double(CK_all{e1,f1});
    end
end

for g1=2:1:n
    for h1=2:1:3
        Cr(g1-1,h1)=str2double(Cr_all{g1,h1});
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
in=MixNumber*m;
a=n;
up=30000000;
low1=[0,0,0,0,0];
low2=[0,0,0,0];
low3=[0,0,0];
upper1=[up,1,up,1,1];
upper2=[up,1,up,1];
upper3=[up,up,1,1];
upper4=[up,up,1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
o = statset('mlecustom');
o.FunValCheck = 'off';
o.MaxIter=1000000;
o.TolFun=0.0001;
simu_data1=zeros(n,m); %simu_data1:K:r=1:9
Kr_simu=zeros(a,b);
for r=1:1:b
    for q=1:1:in
        for p=1:1:n
            judge=isnan(CK(p,3));
            if judge==0
                simu_data1(p,q)=nbinrnd(CK(p,2),CK(p,3));
            else
                simu_data1(p,q)=poissrnd(CK(p,2));
            end
        end
    end
    for q=in+2:1:m+1
        for p=1:1:n
            
            judge=isnan(Cr(p,3));
            if judge==0
                simu_data1(p,q)=nbinrnd(Cr(p,2),Cr(p,3));
            else
                simu_data1(p,q)=poissrnd(Cr(p,2));
            end
        end
    end
    for s=1:1:a
        x=simu_data1(s,:);
        x1=x(2:end);
        xK=x1(1:(length(x1)*(in/m)));
        xr=x1(((length(x1)*(in/m))+1):length(x1));
        variation_K=var(xK);
        average_K=mean(xK);
        variation_r=var(xr);
        average_r=mean(xr);
        vK=roundn(variation_K,-10);
        aK=roundn(average_K,-10);
        vr=roundn(variation_r,-10);
        ar=roundn(average_r,-10);
        if vK > aK
            if vr > ar
                parmhat_K=nbinfit(xK);
                parmhat_r=nbinfit(xr);
                mixedpdf=@(x,r1,p1,r2,p2,rho)(rho*nbinpdf(x,r1,p1)+(1-rho)*nbinpdf(x,r2,p2));
                phat1=mle(x1,'pdf',mixedpdf,'start',[parmhat_K(1),parmhat_K(2),parmhat_r(1),parmhat_r(2),MixNumber],'options',o,'lowerbound',low1,'upperbound',upper1);
                phat1=mle(x1,'pdf',mixedpdf,'start',phat1,'options',o,'lowerbound',low1,'upperbound',upper1);
                phat1=mle(x1,'pdf',mixedpdf,'start',phat1,'options',o,'lowerbound',low1,'upperbound',upper1);
                mf=@(x)(phat1(5)*nbinpdf(x,phat1(1),phat1(2))+(1-phat1(5))*nbinpdf(x,phat1(3),phat1(4)));
                randm=rand;
                cdfvalue=mf(0);
                if cdfvalue>randm
                    Kr_simu(s,r)=0;
                else
                    for x=1:1:10000000
                        cdfvalue=cdfvalue+mf(x);
                        if cdfvalue>randm
                            Kr_simu(s,r)=x;
                            break
                        end
                    end
                end
            elseif vr < ar
                parmhat_K=nbinfit(xK);
                lambdahat_r=poissfit(xr);
                mixedpdf=@(x,r,p,lambda,rho)(rho*nbinpdf(x,r,p)+(1-rho)*poisspdf(x,lambda));
                phat1=mle(x1,'pdf',mixedpdf,'start',[parmhat_K(1),parmhat_K(2),lambdahat_r,MixNumber],'options',o,'lowerbound',low2,'upperbound',upper2);
                phat1=mle(x1,'pdf',mixedpdf,'start',phat1,'options',o,'lowerbound',low2,'upperbound',upper2);
                phat1=mle(x1,'pdf',mixedpdf,'start',phat1,'options',o,'lowerbound',low2,'upperbound',upper2);
                mf=@(x)(phat1(4)*nbinpdf(x,phat1(1),phat1(2))+(1-phat1(4))*poisspdf(x,phat1(3)));
                randm=rand;
                cdfvalue=mf(0);
                if cdfvalue>randm
                    Kr_simu(s,r)=0;
                else
                    for x=1:1:10000000
                        cdfvalue=cdfvalue+mf(x);
                        if cdfvalue>randm
                            Kr_simu(s,r)=x;
                            break
                        end
                    end
                end
            elseif vr == ar
                parmhat_K=nbinfit(xK);
                lambdahat_r=poissfit(xr);
                mixedpdf=@(x,r,p,lambda,rho)(rho*nbinpdf(x,r,p)+(1-rho)*poisspdf(x,lambda));
                phat1=mle(x1,'pdf',mixedpdf,'start',[parmhat_K(1),parmhat_K(2),lambdahat_r,MixNumber],'options',o,'lowerbound',low2,'upperbound',upper2);
                phat1=mle(x1,'pdf',mixedpdf,'start',phat1,'options',o,'lowerbound',low2,'upperbound',upper2);
                phat1=mle(x1,'pdf',mixedpdf,'start',phat1,'options',o,'lowerbound',low2,'upperbound',upper2);
                mf=@(x)(phat1(4)*nbinpdf(x,phat1(1),phat1(2))+(1-phat1(4))*poisspdf(x,phat1(3)));
                randm=rand;
                cdfvalue=mf(0);
                if cdfvalue>randm
                    Kr_simu(s,r)=0;
                else
                    for x=1:1:10000000
                        cdfvalue=cdfvalue+mf(x);
                        if cdfvalue>randm
                            Kr_simu(s,r)=x;
                            break
                        end
                    end
                end
            end
        elseif vK < aK
            if vr < ar
                lambdahat_K=poissfit(xK);
                lambdahat_r=poissfit(xr);
                mixedpdf=@(x,lambda1,lambda2,rho)(rho*poisspdf(x,lambda1)+(1-rho)*poisspdf(x,lambda2));
                phat1=mle(x1,'pdf',mixedpdf,'start',[lambdahat_K,lambdahat_r,MixNumber],'options',o,'lowerbound',low3,'upperbound',upper4);
                phat1=mle(x1,'pdf',mixedpdf,'start',phat1,'options',o,'lowerbound',low3,'upperbound',upper4);
                phat1=mle(x1,'pdf',mixedpdf,'start',phat1,'options',o,'lowerbound',low3,'upperbound',upper4);
                mf=@(x)(phat1(3)*poisspdf(x,phat1(1))+(1-phat1(3))*poisspdf(x,phat1(2)));
                randm=rand;
                cdfvalue=mf(0);
                if cdfvalue>randm
                    Kr_simu(s,r)=0;
                else
                    for x=1:1:10000000
                        cdfvalue=cdfvalue+mf(x);
                        if cdfvalue>randm
                            Kr_simu(s,r)=x;
                            break
                        end
                    end
                end
            elseif vr > ar
                lambdahat_K=poissfit(xK);
                parmhat_r=nbinfit(xr);
                mixedpdf=@(x,lambda,r,p,rho)(rho*poisspdf(x,lambda)+(1-rho)*nbinpdf(x,r,p));
                phat1=mle(x1,'pdf',mixedpdf,'start',[lambdahat_K,parmhat_r(1),parmhat_r(2),MixNumber],'options',o,'lowerbound',low2,'upperbound',upper3);
                phat1=mle(x1,'pdf',mixedpdf,'start',phat1,'options',o,'lowerbound',low2,'upperbound',upper3);
                phat1=mle(x1,'pdf',mixedpdf,'start',phat1,'options',o,'lowerbound',low2,'upperbound',upper3);
                mf=@(x)(phat1(4)*poisspdf(x,phat1(1))+(1-phat1(4))*nbinpdf(x,phat1(2),phat1(3)));
                randm=rand;
                cdfvalue=mf(0);
                if cdfvalue>randm
                    Kr_simu(s,r)=0;
                else
                    for x=1:1:10000000
                        cdfvalue=cdfvalue+mf(x);
                        if cdfvalue>randm
                            Kr_simu(s,r)=x;
                            break
                        end
                    end
                end
            elseif vr == ar
                lambdahat_K=poissfit(xK);
                lambdahat_r=poissfit(xr);
                mixedpdf=@(x,lambda1,lambda2,rho)(rho*poisspdf(x,lambda1)+(1-rho)*poisspdf(x,lambda2));
                phat1=mle(x1,'pdf',mixedpdf,'start',[lambdahat_K,lambdahat_r,MixNumber],'options',o,'lowerbound',low3,'upperbound',upper4);
                phat1=mle(x1,'pdf',mixedpdf,'start',phat1,'options',o,'lowerbound',low3,'upperbound',upper4);
                phat1=mle(x1,'pdf',mixedpdf,'start',phat1,'options',o,'lowerbound',low3,'upperbound',upper4);
                mf=@(x)(phat1(3)*poisspdf(x,phat1(1))+(1-phat1(3))*poisspdf(x,phat1(2)));
                randm=rand;
                cdfvalue=mf(0);
                if cdfvalue>randm
                    Kr_simu(s,r)=0;
                else
                    for x=1:1:10000000
                        cdfvalue=cdfvalue+mf(x);
                        if cdfvalue>randm
                            Kr_simu(s,r)=x;
                            break
                        end
                    end
                end
            end
        elseif vK == aK
            if vr < ar
                lambdahat_K=poissfit(xK);
                lambdahat_r=poissfit(xr);
                mixedpdf=@(x,lambda1,lambda2,rho)(rho*poisspdf(x,lambda1)+(1-rho)*poisspdf(x,lambda2));
                phat1=mle(x1,'pdf',mixedpdf,'start',[lambdahat_K,lambdahat_r,MixNumber],'options',o,'lowerbound',low3,'upperbound',upper4);
                phat1=mle(x1,'pdf',mixedpdf,'start',phat1,'options',o,'lowerbound',low3,'upperbound',upper4);
                phat1=mle(x1,'pdf',mixedpdf,'start',phat1,'options',o,'lowerbound',low3,'upperbound',upper4);
                mf=@(x)(phat1(3)*poisspdf(x,phat1(1))+(1-phat1(3))*poisspdf(x,phat1(2)));
                randm=rand;
                cdfvalue=mf(0);
                if cdfvalue>randm
                    Kr_simu(s,r)=0;
                else
                    for x=1:1:10000000
                        cdfvalue=cdfvalue+mf(x);
                        if cdfvalue>randm
                            Kr_simu(s,r)=x;
                            break
                        end
                    end
                end
            elseif vr > ar
                lambdahat_K=poissfit(xK);
                parmhat_r=nbinfit(xr);
                mixedpdf=@(x,lambda,r,p,rho)(rho*poisspdf(x,lambda)+(1-rho)*nbinpdf(x,r,p));
                phat1=mle(x1,'pdf',mixedpdf,'start',[lambdahat_K,parmhat_r(1),parmhat_r(2),MixNumber],'options',o,'lowerbound',low2,'upperbound',upper3);
                phat1=mle(x1,'pdf',mixedpdf,'start',phat1,'options',o,'lowerbound',low2,'upperbound',upper3);
                phat1=mle(x1,'pdf',mixedpdf,'start',phat1,'options',o,'lowerbound',low2,'upperbound',upper3);
                mf=@(x)(phat1(4)*poisspdf(x,phat1(1))+(1-phat1(4))*nbinpdf(x,phat1(2),phat1(3)));
                randm=rand;
                cdfvalue=mf(0);
                if cdfvalue>randm
                    Kr_simu(s,r)=0;
                else
                    for x=1:1:10000000
                        cdfvalue=cdfvalue+mf(x);
                        if cdfvalue>randm
                            Kr_simu(s,r)=x;
                            break
                        end
                    end
                end
            elseif vr == ar
                lambdahat_K=poissfit(xK);
                lambdahat_r=poissfit(xr);
                mixedpdf=@(x,lambda1,lambda2,rho)(rho*poisspdf(x,lambda1)+(1-rho)*poisspdf(x,lambda2));
                phat1=mle(x1,'pdf',mixedpdf,'start',[lambdahat_K,lambdahat_r,MixNumber],'options',o,'lowerbound',low3,'upperbound',upper4);
                phat1=mle(x1,'pdf',mixedpdf,'start',phat1,'options',o,'lowerbound',low3,'upperbound',upper4);
                phat1=mle(x1,'pdf',mixedpdf,'start',phat1,'options',o,'lowerbound',low3,'upperbound',upper4);
                mf=@(x)(phat1(3)*poisspdf(x,phat1(1))+(1-phat1(3))*poisspdf(x,phat1(2)));
                randm=rand;
                cdfvalue=mf(0);
                if cdfvalue>randm
                    Kr_simu(s,r)=0;
                else
                    for x=1:1:10000000
                        cdfvalue=cdfvalue+mf(x);
                        if cdfvalue>randm
                            Kr_simu(s,r)=x;
                            break
                        end
                    end
                end
            end
        end
        r
        s
    end
end
clear simu_data1;
Kr_simu_results=cell(n+1,b+1);
for j=1:1:b+1
    for i=1:1:n+1
        if j==1
            Kr_simu_results{i,j}=CK_all{i,j};
        else
            if i==1
                Kr_simu_results{i,j}=(dirname);
            else
                Kr_simu_results{i,j}=Kr_simu(i-1,j-1);
            end
        end
    end
end
clear Kr_simu;
save (['./' dirname '/Kr_simu_results'], 'Kr_simu_results');
fid1=fopen(['./' dirname '/Kr_simu_results.txt'],'wt');
for r=1
    for s=1:1:b+1
        if s<b+1
            fprintf(fid1,'%s\t',Kr_simu_results{r,s});
        else
            fprintf(fid1,'%s\n',Kr_simu_results{r,s});
        end
    end
end
for r=2:1:n+1
    for s=1:1:b+1     
            if s==1
                fprintf(fid1,'%s\t',Kr_simu_results{r,s});
            elseif (s==b+1)
                fprintf(fid1,'%g\n',Kr_simu_results{r,s});
            else
                fprintf(fid1,'%g\t',Kr_simu_results{r,s});
            end
    end
end
fclose(fid1);
clear Kr_simu_results;
end