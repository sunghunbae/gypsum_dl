# Copyright 2023 Jacob D. Durrant
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


# gypsum_dl/gypsum_dl/
# Including the below allows other programs to import functions from
# gypsum-DL.
# import sys
# import os

# current_dir = os.path.dirname(os.path.realpath(__file__))
# gypsum_gypsum_dir = current_dir
# gypsum_top_dir = os.path.dirname(gypsum_gypsum_dir)
# sys.path.extend([gypsum_gypsum_dir, gypsum_top_dir])

__version__ = '1.2.1'

from .molstates import MolStates