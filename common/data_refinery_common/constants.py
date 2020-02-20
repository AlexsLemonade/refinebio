from data_refinery_common.utils import get_env_variable

LOCAL_ROOT_DIR = get_env_variable("LOCAL_ROOT_DIR", "/home/user/data_store")
# We store what salmon ouptuts as its version, therefore for
# comparisions or defaults we shouldn't just store the version string,
# we need something with the pattern: 'salmon X.X.X'
CURRENT_SALMON_VERSION = "salmon " + get_env_variable("SALMON_VERSION", "0.13.1")
CHUNK_SIZE = 1024 * 256  # chunk_size is in bytes
# Let this fail if SYSTEM_VERSION is unset.
SYSTEM_VERSION = get_env_variable("SYSTEM_VERSION")
