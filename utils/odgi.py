from utils.common import run


def ogBuild(gfaFile, ogFile, threads):
    ODGIBuildCmd = (
        ["odgi", "build", "-g", gfaFile, "-s", "-o", ogFile]
        if threads == 1
        else ["odgi", "build", "-g", gfaFile, "-s", "-o", ogFile, "-t", threads]
    )
    run(ODGIBuildCmd)


def ogSort(ogFile, ogSortedFile, threads):
    ODGIBuildCmd = (
        ["odgi", "sort", "-i", ogFile, "-o", ogSortedFile]
        if threads == 1
        else ["odgi", "sort", "-i", ogFile, "-o", ogSortedFile, "-t", threads]
    )
    run(ODGIBuildCmd)


def ogExtract(ogFile, extractogFile, tPath, tRange, threads):
    ODGIExtractCmd = (
        [
            "odgi",
            "extract",
            "-i",
            ogFile,
            "-E",
            "-r",
            tPath + ":" + str(tRange[0]) + "-" + str(tRange[1]),
            "-d",
            "0",
            "-o",
            extractogFile,
        ]
        if threads == 1
        else [
            "odgi",
            "extract",
            "-i",
            ogFile,
            "-E",
            "-r",
            tPath + ":" + str(tRange[0]) + "-" + str(tRange[1]),
            "-d",
            "0",
            "-o",
            extractogFile,
            "-t",
            threads,
        ]
    )
    run(ODGIExtractCmd)


def ogExtractBed(ogFile, extractOgFile, bedFile, threads):
    ODGIExtractCmd = (
        [
            "odgi",
            "extract",
            "-i",
            ogFile,
            "-b",
            bedFile,
            "-d",
            "0",
            "-o",
            extractOgFile,
        ]
        if threads == 1
        else [
            "odgi",
            "extract",
            "-i",
            ogFile,
            "-b",
            bedFile,
            "-d",
            "0",
            "-o",
            extractOgFile,
            "-t",
            threads,
        ]
    )
    run(ODGIExtractCmd)


def ogPath(ogFile, csvFile, threads):
    ODGIPathCmd = (
        ["odgi", "paths", "-i", ogFile, "-H"]
        if threads == 1
        else ["odgi", "paths", "-i", ogFile, "-H"]
    )
    if_success, stdout, stderr = run(ODGIPathCmd)
    with open(csvFile, "w") as f:
        f.write(stdout)


def ogView(ogFile, gfaFile, threads):
    ODGIViewCmd = (
        ["odgi", "view", "-i", ogFile, "-g"]
        if threads == 1
        else ["odgi", "view", "-i", ogFile, "-g", "-t", threads]
    )
    if_success, stdout, stderr = run(ODGIViewCmd)
    with open(gfaFile, "w") as f:
        f.write(stdout)
