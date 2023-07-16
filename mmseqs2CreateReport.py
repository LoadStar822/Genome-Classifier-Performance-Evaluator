import os
import subprocess

directory = "/home/zqtianqinzhong/software/ART/datasets/mmseqs2_results"

for filename in os.listdir(directory):
    if filename.endswith(".out.0"):
        new_filename = filename.rsplit(".", 1)[0]

        command = "mmseqs taxonomyreport /home/zqtianqinzhong/software/mmseqs2/nr {resultDB} {taxonomyResult_report}"

        filled_command = command.format(resultDB=os.path.join(directory, new_filename),
                                        taxonomyResult_report=os.path.join(directory, new_filename))

        # 执行命令
        subprocess.run(filled_command, shell=True)
