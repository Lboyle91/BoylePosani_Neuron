function [full,vidframes_bytrial,vidframesall]=CorrectDroppedFrames(initial, dropped, numtests, vidframesall,vidframes_bytrial)

last=0;
vf_bt_static=vidframes_bytrial;


full=[];
    for i=1:numtests
        first=last+1;
        currvidfr=vidframes_bytrial{i};
        currdroppedtest=dropped{i};
        
        currvidfrstat=vf_bt_static{i};
        ct=initial(:,currvidfrstat(1):currvidfrstat(2));
        currtrifr=vidframesall{i};
        currtrifrstatic=currtrifr;
        
        for test=1:length(currdroppedtest)
            currfr=currtrifr{test};
            currfrstat=currtrifrstatic{test};
            
            first=currfr(1);
        	last=currfr(2);
            
            tr1=currfrstat(1);
            tr2=currfrstat(2);
            currtrace=ct(:,tr1:tr2);
            currdropped=currdroppedtest{test};
            if ~isnan(currdropped(1))
                nwtrace=currtrace;
                for dropidx=1:length(currdropped)
                    traceA=nwtrace(:,1:currdropped(dropidx)-1);
                    traceB=nwtrace(:,currdropped(dropidx):length(nwtrace));
                    nwtrace=[traceA NaN(length(initial(:,1)),1) traceB];
                end
                last=last+length(currdropped); %add change to current frame indices
                currtrifr{test}=[first last];
                
                for a=(test+1):length(currtrifr) %add change to indices within same test
                    currvfa=currtrifr{a};
                    currvfa=currvfa+length(currdropped);
                    currtrifr{a}=currvfa;
                end
                vidframesall{i}=currtrifr; 
                
                currvf=vidframes_bytrial{i};
                vidframes_bytrial{i}=[currvf(1) currvf(2)+length(currdropped)]; %add change to current global test indices
                
                for b=(i+1):numtests %add change to future test indices
                    cur=vidframes_bytrial{b}; %add change to global ind
                    cur=cur+length(currdropped);
                    vidframes_bytrial{b}=cur;
                end
            else
                nwtrace=currtrace;
            end
            full=[full nwtrace];
        end
    end

        