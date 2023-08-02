start_logging() {

    ACTION=$1
    START_TIME_HHMM=$(date +" %T")
    START_TIME=$(date +%s)
}

stop_logging() {

    END_TIME_HHMM=$(date +"%T")
    END_TIME=$(date +%s)
    ELAPSED=$(expr $END_TIME - $START_TIME)
    DURATION=$(date -u -d @${ELAPSED} +"%T")
    printf '%-80s   %8s   %8s   %8s\n' "${ACTION}" ${START_TIME_HHMM} ${END_TIME_HHMM} ${DURATION} >>${LOG_FILE}
}

