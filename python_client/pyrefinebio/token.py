import os

import yaml
from pyrefinebio.http import get, post, put


class Token:
    """
        Handles the creation, activation, saving, and loading of api tokens.

        These tokens can be used in requests that provide urls to download computed files.

        Please review refine.bio's [Terms of Use](https://www.refine.bio/terms) and [Privacy Policy](https://www.refine.bio/privacy)
        before use of these tokens.
    """

    def create_token(self, email_address):
        """
            creates a token and emails the terms and conditions to a specified email.

            parameters:

                email_address (str): the email that the terms and conditions should be sent to.
        """

        response = post("token")

        token_id = response["id"]

        if email_address:
            terms = response["terms_and_conditions"]
            # email terms

        return token_id

    def agree_to_terms_and_conditions(self, api_token):
        """
            Activates a token.

            Activating a token indicates agreement with refine.bio's
            [Terms of Use](https://www.refine.bio/terms) and
            [Privacy Policy](https://www.refine.bio/privacy).

            parameters:

                api_token (str): the uuid string identifying the token
                                 you want to activate.
        """

        response = put("token" + api_token, payload={"is_activated": True})
        return response

    def save_token(self, api_token, file_path=os.getenv("CONFIG_FILE", "~/.refinebio.yaml")):
        """
            Saves a token to a file.

            parameters:

                api_token (str): the uuid string identifying the token
                                 you want to save.
                file_path (str): the path to the file where the token should be
                                 saved. Defaults to `~/.refinebio.yaml`. Alternatively
                                 you can set this path in an environment variable `CONFIG_FILE`.
        """

        with open(file_path, "w+") as file:
            yaml.dump({"token": api_token}, file)

    def load_token(self, file_path=os.getenv("CONFIG_FILE", "~/.refinebio.yaml")):
        """
            Loads a token from a file.

            parameters:

                file_path (str): the path to the file where the token should be
                                 loaded from. Defaults to `~/.refinebio.yaml`. Alternatively
                                 you can set this path in an environment variable `CONFIG_FILE`.
        """

        with open(file_path) as file:
            tokens = yaml.full_load(file)
            return tokens["token"]
