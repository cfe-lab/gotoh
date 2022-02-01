""" Define the configuration parameters of the alignment webserver
"""


# from enum import Enum
from pydantic import BaseModel


class WebServerConfig(BaseModel):
    """Define the configuration parameters of the webserver."""

    # parameters for creating API client tokens
    # to get a string like this run:
    # openssl rand -hex 32
    SECRET_KEY: str = "my_secret_key"
    HASH_ALGORITHM: str = "HS256"
    ACCESS_TOKEN_EXPIRE_MINUTES: int = 15

    # if provided, this is a YAML file from which
    # a reportparams.GlobalConfig instance will be retrieved.
    # Leave this black to produce a default GlobalConfig istance
    # using self.get_global_config() below..
    global_config_file: str = None

    # email notification configuration
    SMTP_SERVER: str = ''
    SMTP_PORT: int = 25
    # This person will be sent an email when an automatically started
    # monthly report is successfully completed
    AUTO_EMAIL_USER: str = ''
    # This person will be sent an email whenever an error requiring
    # developer attention occurs.
    DEVEL_EMAIL_USER: str = ''
