# hindcast track paths -- source AFTER env.sh
if [ -z "${CROCO_ROOT}" ]; then
    echo "!! CROCO_ROOT is not set -- run 'source ~/croco_work/env.sh' first"
    return 1 2>/dev/null || exit 1
fi
export CROCO_CONFIGS_ROOT=${CROCO_ROOT}/hindcast/configs
export CROCO_RUNS_ROOT=${CROCO_ROOT}/hindcast/scratch
echo "  track = hindcast  (configs: hindcast/configs, runs: hindcast/scratch)"
