/**
 *  Copyright 2007-2008 University Of Southern California
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *  http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing,
 *  software distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */

import java.io.BufferedReader;
//import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.InputStreamReader;
import java.io.IOException;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.StringTokenizer;
import edu.isi.pegasus.planner.dax.*;
import edu.isi.pegasus.planner.dax.Invoke.WHEN;

public class cwb_dag2dax {

  /**
   * Create CWB DAX
   * @param args
   */

  public static void main(String[] args) {
    if (args.length != 6) {
      String usage = "Usage: ";
      usage += "java ADAG <pegasus_bin_dir> <filename.dag> <filename.tgz> ";
      usage += "<filename.dax> <filename.out> <filename.watenv>";
      System.out.println(usage);
      System.exit(1);
    }

    try {
      CWB(args[0], args[1], args[2], args[4], args[5]).writeToFile(args[3]);
    }
    catch (Exception e) {
      e.printStackTrace();
    }
  }

  private static String getValue(String s) {
    Pattern p = Pattern.compile("\"([^\"]*)\"");
    Matcher m = p.matcher(s);
    String result = null;
    if (m.find()) {
      result = m.group(1);
    }
    return result;
  } 

  private static ADAG CWB(String pegasus_bin_dir, String dagFile, String tgzFile, 
                          String outFile, String watenvFile) throws Exception {

    java.io.File cwdFile = new java.io.File ("..");
    String cwd = cwdFile.getCanonicalPath();
    String output_dir = "output"; // warning : this name is defined in the $CWB_PARAMETERS_FILE
    String condor_dir = "condor"; // warning : this name is defined in the $CWB_PARAMETERS_FILE
    String cwb_exec   = "cwb.sh"; // warning : this name is defined in cwb_mkdir.C

    ADAG dax = new ADAG("cwb");

    // make full tgzFile path
    File workingdir = new File(tgzFile);
    workingdir.addPhysicalFile("file://"+cwd+"/"+condor_dir+"/"+tgzFile, "local");
    dax.addFile(workingdir);

    // make full cwb_exec path
    Executable cwb_net = new Executable("cwb_net", "cwb_net", "4.0");
    cwb_net.setArchitecture(Executable.ARCH.X86).setOS(Executable.OS.LINUX);
    cwb_net.setInstalled(false);
    cwb_net.addPhysicalFile("file://"+cwd+"/"+condor_dir+"/"+cwb_exec, "local");

    dax.addExecutable(cwb_net);

    // Open dagFile 
    FileInputStream fstream = new FileInputStream(dagFile);
    BufferedReader br = new BufferedReader(new InputStreamReader(fstream));
 
    String strLine = null;
    String jobName = null;
    String jobId = null;
    String cwbUFile = null;
    String cwbStage = null;
    StringTokenizer st = null;
    Job job = null;
    File output = null;
    Pattern p = Pattern.compile("\"([^\"]*)\"");
    Matcher m = null;
 
    //Read dagFile Line By Line
    while ((strLine = br.readLine()) != null)   {
      if (strLine.startsWith("VARS")) {
        st = new StringTokenizer(strLine.trim(), " ");
        st.nextToken();
        jobName = st.nextToken();        
        jobId = getValue(st.nextToken());
        cwbUFile = getValue(st.nextToken());
        cwbStage = getValue(st.nextToken());
        output = new File(outFile + ".job" + jobId + ".tgz");
        output.setRegister(false);

        job = new Job(jobName, "cwb_net", "cwb_net", "4.0");
        job.addArgument(jobId);
        job.addArgument(cwbUFile);
        job.addArgument(cwbStage);
        job.addArgument(cwdFile.getCanonicalFile().getName());
        job.addArgument(tgzFile);
        job.addArgument(outFile);
        job.addArgument(watenvFile);
        job.setStdout(outFile + ".job" + jobId + ".txt", File.TRANSFER.OPTIONAL, false);
        //job.setStderr("stderr_stage1.txt", File.TRANSFER.OPTIONAL, false);
        job.uses(workingdir, File.LINK.INPUT, File.TRANSFER.TRUE, false);
        job.uses(output, File.LINK.OUTPUT, File.TRANSFER.OPTIONAL, false);

        String arguments = "VARS A" + jobId + " PID=\\\"" + jobId + "\\\" CWB_UFILE=\\\"" + 
                           cwbUFile + "\\\" CWB_STAGE=\\\"" + cwbStage + "\\\"\n";

        // multistage analysis : add root job file 
        if (cwbUFile.endsWith(".root")) {
          File uFile = null;

          if (cwbUFile.startsWith(output_dir+"/")) {
            uFile = new File(cwbUFile.substring(output_dir.length()+1));
          } else {
            uFile = new File(cwbUFile);
          }

          uFile.addPhysicalFile("file://"+cwd+"/"+cwbUFile, "local");

          job.uses(uFile, File.LINK.INPUT, File.TRANSFER.TRUE, false);
          dax.addFile(uFile);
        }

        dax.addJob(job);

        // at the end of job uncompress tgz in the local working directory
        String tar_cmd = "/bin/tar -xzvf "+cwd+"/"+output_dir+"/"+output.getName()+" -C "+cwd;
        dax.addInvoke(WHEN.at_end, tar_cmd);

        // at the end of job remove tgz 
        String rm_cmd = "/bin/rm -f "+cwd+"/"+output_dir+"/"+output.getName();
        dax.addInvoke(WHEN.at_end, rm_cmd);

        // at the end move the .txt files from the output dir to the log dir
        String mv_cmd = "/bin/mv "+cwd+"/"+output_dir+"/"+outFile+".job"+jobId+".txt "+cwd+"/log/";
        dax.addInvoke(WHEN.at_end, mv_cmd);
      }
    }

    //Close the input stream
    br.close();

    return dax;
  }
}
