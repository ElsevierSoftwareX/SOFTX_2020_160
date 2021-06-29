#!/bin/tcsh -f

condor_q -submitter $USER -g \
-format "%d\t" JobStatus \
-format "%s\t" DAGNodeName \
-format "%d\t" EnteredCurrentStatus \
-format "%d\t" ServerTime \
-format "%d\t" ClusterId \
-format "%d\t" ImageSize \
-format "%s\n" Cmd | awk '
BEGIN {
        nR=0;
        nH=0;
        nI=0;
        nJ=0;

        USER = ENVIRON["USER"];
        PWD = ENVIRON["PWD"];
	PWD=sprintf("%s%s",PWD,"/");

        printf("----------------------------------------------------------------------------\n");
        printf("User : %s\n",USER);
        #printf("Working Dir : %s\n",PWD);
        printf("----------------------------------------------------------------------------\n");
	printf("JOB-ID\tJOB-TYPE\tJOB-STATUS\tRUN-TIME\tMEM(KB)\tCLUSTER-ID\n"); 
        printf("---------------------------------------------------------------------------\n");
}
{
	Cmd = $7; 
        where = match(Cmd, PWD);
        if(where > 0) {

		DAGNodeName = $2; 
		EnteredCurrentStatus = $3; 
		ServerTime = $4; 
		ClusterId = $5; 
		ImageSize = $6; 

                if(match(Cmd, "cwb.sh")>0) JobType="net";
                if(match(Cmd, "ced.sh")>0) JobType="ced";

                nJ++;

		if($1 == 1) {JobStatus="Idle";nI++;} 
		if($1 == 2) {JobStatus="Running";nR++} 
		if($1 == 5) {JobStatus="Hold";nH++;} 
		RunTime = strftime("%H:%M:%S", ServerTime - EnteredCurrentStatus); 
		RunDays = (ServerTime - EnteredCurrentStatus) / 86400;   

		DAGNodeName=sprintf("%s",substr(DAGNodeName,2)); 

		printf("%s\t%s\t\t%s\t\t%1d+%s\t%d\t%s\n",DAGNodeName,JobType,JobStatus,RunDays,RunTime,ImageSize/1024.,ClusterId); 
        }
}
END {
        printf("----------------------------------------------------------------------------\n");
        printf("%d jobs; %d idle, %d running, %d held\n",nJ,nI,nR,nH);
        printf("----------------------------------------------------------------------------\n");
}
'
