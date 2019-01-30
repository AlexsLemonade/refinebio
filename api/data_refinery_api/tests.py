import json
import random
import time

from django.contrib.auth.models import User
from django.http import HttpResponseForbidden, HttpResponseServerError
from django.urls import reverse
from rest_framework import status
from rest_framework.test import APITestCase
from unittest.mock import patch

from data_refinery_api.serializers import (
    DetailedExperimentSerializer,
    DetailedSampleSerializer,
    ExperimentSerializer,
    InstitutionSerializer,
    OrganismIndexSerializer,
    OrganismSerializer,
    PlatformSerializer,
    SampleSerializer,

    # Jobs
    DownloaderJobSerializer,
    ProcessorJobSerializer,
    ProcessorSerializer,
    SurveyJobSerializer,
)
from data_refinery_api.views import ExperimentList
from data_refinery_common.utils import get_env_variable
from data_refinery_common.models import (
    ComputationalResult,
    ComputedFile,
    Dataset,
    DownloaderJob,
    DownloaderJobOriginalFileAssociation,
    Experiment,
    ExperimentAnnotation,
    ExperimentOrganismAssociation,
    ExperimentSampleAssociation,
    Organism,
    OrganismIndex,
    OriginalFile,
    OriginalFileSampleAssociation,
    Processor,
    ProcessorJob,
    ProcessorJobOriginalFileAssociation,
    Sample,
    SampleAnnotation,
    SampleResultAssociation,
)
from data_refinery_common.models.documents import (
    ExperimentDocument
)

class APITestCases(APITestCase):
    def setUp(self):
        # Saving this for if we have protected endpoints
        # self.superuser = User.objects.create_superuser('john', 'john@snow.com', 'johnpassword')
        # self.client.login(username='john', password='johnpassword')
        # self.user = User.objects.create(username="mike")

        experiment = Experiment()
        experiment.accession_code = "GSE000"
        experiment.title = "NONONONO"
        experiment.description = "Boooooourns. Wasabi."
        experiment.technology = "RNA-SEQ"
        experiment.save()

        experiment = Experiment()
        experiment.accession_code = "GSE123"
        experiment.title = "Hey Ho Let's Go"
        experiment.description = "This is a very exciting test experiment. Faygo soda. Blah blah blah."
        experiment.technology = "MICROARRAY"
        experiment.save()
        self.experiment = experiment

        experiment_annotation = ExperimentAnnotation()
        experiment_annotation.data = {"hello": "world", "123": 456}
        experiment_annotation.experiment = experiment
        experiment_annotation.save()

        sample = Sample()
        sample.title = "123"
        sample.accession_code = "123"
        sample.is_processed = True
        sample.save()

        sample = Sample()
        sample.title = "789"
        sample.accession_code = "789"
        sample.is_processed = True
        sample.save()
        self.sample = sample

        sample_annotation = SampleAnnotation()
        sample_annotation.data = {"goodbye": "world", "789": 123}
        sample_annotation.sample = sample
        sample_annotation.save()

        original_file = OriginalFile()
        original_file.save()

        original_file_sample_association = OriginalFileSampleAssociation()
        original_file_sample_association.sample = sample
        original_file_sample_association.original_file = original_file
        original_file_sample_association.save()

        downloader_job = DownloaderJob()
        downloader_job.save()

        download_assoc = DownloaderJobOriginalFileAssociation()
        download_assoc.original_file = original_file
        download_assoc.downloader_job = downloader_job
        download_assoc.save()

        processor_job = ProcessorJob()
        processor_job.save()

        processor_assoc = ProcessorJobOriginalFileAssociation()
        processor_assoc.original_file = original_file
        processor_assoc.processor_job = processor_job
        processor_assoc.save()

        experiment_sample_association = ExperimentSampleAssociation()
        experiment_sample_association.sample = sample
        experiment_sample_association.experiment = experiment
        experiment_sample_association.save()

        result = ComputationalResult()
        result.save()

        sra = SampleResultAssociation()
        sra.sample = sample
        sra.result = result
        sra.save()

        result = ComputationalResult()
        result.save()

        sra = SampleResultAssociation()
        sra.sample = sample
        sra.result = result
        sra.save()

        return

    def tearDown(self):
        """ Good bye """
        Experiment.objects.all().delete()
        ExperimentAnnotation.objects.all().delete()
        Sample.objects.all().delete()
        SampleAnnotation.objects.all().delete()
        Sample.objects.all().delete()
        SampleAnnotation.objects.all().delete()

    def test_all_endpoints(self):
        response = self.client.get(reverse('experiments'))
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(response['X-Source-Revision'], get_env_variable('SYSTEM_VERSION'))

        response = self.client.get(reverse('experiments'), kwargs={'page': 1})
        self.assertEqual(response.status_code, status.HTTP_200_OK)

        response = self.client.get(reverse('experiments_detail', kwargs={'pk': '1'}))
        self.assertEqual(response.status_code, status.HTTP_200_OK)

        response = self.client.get(reverse('samples'))
        self.assertEqual(response.status_code, status.HTTP_200_OK)

        response = self.client.get(reverse('samples'), {'ids': str(self.sample.id) + ',1000'})
        self.assertEqual(response.status_code, status.HTTP_200_OK)

        response = self.client.get(reverse('samples'), {'accession_codes': str(self.sample.accession_code) + ',1000'})
        self.assertEqual(response.status_code, status.HTTP_200_OK)

        response = self.client.get(reverse('samples'), kwargs={'page': 1})
        self.assertEqual(response.status_code, status.HTTP_200_OK)

        response = self.client.get(reverse('samples_detail', kwargs={'pk': '1'}))
        self.assertEqual(response.status_code, status.HTTP_200_OK)

        response = self.client.get(reverse('organisms'))
        self.assertEqual(response.status_code, status.HTTP_200_OK)

        response = self.client.get(reverse('platforms'))
        self.assertEqual(response.status_code, status.HTTP_200_OK)

        response = self.client.get(reverse('institutions'))
        self.assertEqual(response.status_code, status.HTTP_200_OK)

        response = self.client.get(reverse('jobs'))
        self.assertEqual(response.status_code, status.HTTP_200_OK)

        response = self.client.get(reverse('survey_jobs'))
        self.assertEqual(response.status_code, status.HTTP_200_OK)

        response = self.client.get(reverse('downloader_jobs'))
        self.assertEqual(response.status_code, status.HTTP_200_OK)

        response = self.client.get(reverse('processor_jobs'))
        self.assertEqual(response.status_code, status.HTTP_200_OK)

        response = self.client.get(reverse('stats'))
        self.assertEqual(response.status_code, status.HTTP_200_OK)

        response = self.client.get(reverse('results'))
        self.assertEqual(response.status_code, status.HTTP_200_OK)

        response = self.client.get(reverse('results'), kwargs={'page': 1})
        self.assertEqual(response.status_code, status.HTTP_200_OK)

        response = self.client.get(reverse('api_root'))
        self.assertEqual(response.status_code, status.HTTP_200_OK)

        response = self.client.get(reverse('search'))
        self.assertEqual(response.status_code, status.HTTP_200_OK)

        response = self.client.get(reverse('dataset_root'))
        self.assertEqual(response.status_code, status.HTTP_200_OK)

        response = self.client.get(reverse('create_dataset'))
        self.assertEqual(response.status_code, status.HTTP_405_METHOD_NOT_ALLOWED)

        response = self.client.get(reverse('token'))
        self.assertEqual(response.status_code, status.HTTP_200_OK)

    def test_sample_pagination(self):

        response = self.client.get(reverse('samples'))
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(len(response.json()['results']), 2)

        response = self.client.get(reverse('samples'), {'limit': 1})
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(len(response.json()['results']), 1)

        response = self.client.get(reverse('samples'), {'limit': 1, 'order_by': '-title'})
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(response.json()['results'][0]['title'], '789')

        response = self.client.get(reverse('samples'), {'limit': 1, 'order_by': 'title'})
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(response.json()['results'][0]['title'], '123')

    def test_fetching_experiment_samples(self):
        response = self.client.get(reverse('samples'), {'experiment_accession_code': self.experiment.accession_code})
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(len(response.json()['results']), 1)
        self.assertEqual(response.json()['results'][0]['accession_code'], '789')

        # Expect 404 if the experiment accession code isn't valid
        response = self.client.get(reverse('samples'), {'experiment_accession_code': 'wrong-accession-code'})
        self.assertEqual(response.status_code, 404)
        
    def test_compendia(self):
        homo_sapiens = Organism.get_object_for_name("HOMO_SAPIENS")
        danio_rerio = Organism.get_object_for_name("DANIO_RERIO")
        
        result = ComputationalResult()
        result.save()

        hsc1 = ComputedFile()
        hsc1.absolute_file_path = '/null/1.tsv'
        hsc1.filename = '1.tsv'
        hsc1.sha1 = "abc"
        hsc1.size_in_bytes = 1
        hsc1.is_smashable = False
        hsc1.is_qn_target = False
        hsc1.result = result
        hsc1.is_compendia = True
        hsc1.compendia_organism = homo_sapiens
        hsc1.compendia_version = 1
        hsc1.s3_bucket = "dr-compendia"
        hsc1.s3_key = "hsc1.tsv"
        hsc1.save()

        hsc2 = ComputedFile()
        hsc2.absolute_file_path = '/null/2.tsv'
        hsc2.filename = '2.tsv'
        hsc2.sha1 = "abc"
        hsc2.size_in_bytes = 1
        hsc2.is_smashable = False
        hsc2.is_qn_target = False
        hsc2.result = result
        hsc2.is_compendia = True
        hsc2.compendia_organism = homo_sapiens
        hsc2.compendia_version = 2
        hsc2.s3_bucket = "dr-compendia"
        hsc2.s3_key = "hsc2.tsv"
        hsc2.save()

        drc1 = ComputedFile()
        drc1.absolute_file_path = '/null/1.tsv'
        drc1.filename = '1.tsv'
        drc1.sha1 = "abc"
        drc1.size_in_bytes = 1
        drc1.is_smashable = False
        drc1.is_qn_target = False
        drc1.result = result
        drc1.is_compendia = True
        drc1.compendia_organism = danio_rerio
        drc1.compendia_version = 1
        drc1.s3_bucket = "dr-compendia"
        drc1.s3_key = "drc2.tsv"
        drc1.save()

        response = self.client.get(reverse('compendia'))
        self.assertEqual(3, len(response.json()))

    def test_search_and_filter(self):

        sample = Sample()
        sample.accession_code = "XXXXXXXXXXXXXXX"
        sample.is_processed = True
        sample.technology = "RNA-SEQ"
        sample.save()

        # Our Docker image doesn't have the standard dict. >=[
        words = ['the', 'of', 'to', 'and', 'a', 'in', 'is', 'it', 'you', 'that', 'he', 'was', 'for', 'on', 'are', 'with', 'as', 'I', 'his', 'they', 'be', 'at', 'one', 'have', 'this', 'from', 'or', 'had', 'by', 'hot', 'word', 'but', 'what', 'some', 'we', 'can', 'out', 'other', 'were', 'all', 'there', 'when', 'up', 'use', 'your', 'how', 'said', 'an', 'each', 'she', 'which', 'do', 'their', 'time', 'if', 'will', 'way', 'about', 'many', 'then', 'them', 'write', 'would', 'like', 'so', 'these', 'her', 'long', 'make', 'thing', 'see', 'him', 'two', 'has', 'look', 'more', 'day', 'could', 'go', 'come', 'did', 'number', 'sound', 'no', 'most', 'people', 'my', 'over', 'know', 'water', 'than', 'call', 'first', 'who', 'may', 'down', 'side', 'been', 'now', 'find', 'any', 'new', 'work', 'part', 'take', 'get', 'place', 'made', 'live', 'where', 'after', 'back', 'little', 'only', 'round', 'man', 'year', 'came', 'show', 'every', 'good', 'me', 'give', 'our', 'under', 'name', 'very', 'through', 'just', 'form', 'sentence', 'great', 'think', 'say', 'help', 'low', 'line', 'differ', 'turn', 'cause', 'much', 'mean', 'before', 'move', 'right', 'boy', 'old', 'too', 'same', 'tell', 'does', 'set', 'three', 'want', 'air', 'well', 'also', 'play', 'small', 'end', 'put', 'home', 'read', 'hand', 'port', 'large', 'spell', 'add', 'even', 'land', 'here', 'must', 'big', 'high', 'such', 'follow', 'act', 'why', 'ask', 'men', 'change', 'went', 'light', 'kind', 'off', 'need', 'house', 'picture', 'try', 'us', 'again', 'animal', 'point', 'mother', 'world', 'near', 'build', 'self', 'earth', 'father', 'head', 'stand', 'own', 'page', 'should', 'country', 'found', 'answer', 'school', 'grow', 'study', 'still', 'learn', 'plant', 'cover', 'food', 'sun', 'four', 'between', 'state', 'keep', 'eye', 'never', 'last', 'let', 'thought', 'citcitree', 'cross', 'farm', 'hhhh', 'start', 'might', 'stosy', 'saw', 'sar', 'sea', 'draw', 'left', 'late', 'run', "don't", 'while', 'press', 'close', 'night', 'real', 'life', 'few', 'north', 'open', 'seemseegether', 'next', 'white', 'chilchin', 'chiln', 'got', 'walk', 'exampexamase', 'paper', 'group', 'always', 'music', 'those', 'botbotark', 'often', 'letter', 'until', 'mile', 'river', 'car', 'feet', 'care', 'second', 'book', 'carry', 'took', 'science', 'eat', 'room', 'friefr', 'bbban', 'idea', 'fish', 'mountain', 'stop', 'once', 'base', 'hear', 'horse', 'cut', 'sure', 'watch', 'color', 'face', 'wood', 'main', 'enough', 'plain', 'girl', 'usual', 'young', 'ready', 'above', 'ever', 'red', 'list', 'though', 'feel', 'talk', 'bird', 'soon', 'body', 'dog', 'family', 'direct', 'pose', 'leave', 'song', 'measure', 'door', 'product', 'black', 'short', 'numeral', 'class', 'wind', 'question', 'happen', 'complete', 'ship', 'area', 'half', 'rock', 'order', 'fire', 'south', 'problem', 'piece', 'told', 'knew', 'pass', 'since', 'top', 'whole', 'king', 'space', 'heard', 'best', 'hour', 'better', 'true', 'during', 'hundred', 'five', 'rrrember', 'step', 'early', 'hold', 'west', 'groundgroterest', 'reach', 'fast', 'verb', 'sing', 'llsten', 'six', 'table', 'travel', 'less', 'morning', 'ten', 'simple', 'several', 'vowel', 'toward', 'war', 'lay', 'against', 'pattern', 'slow', 'center', 'love', 'person', 'money', 'serve', 'appear', 'road', 'map', 'rain', 'rule', 'govern', 'pull', 'cold', 'notice', 'voice', 'unit', 'powepotown', 'fine', 'certain', 'flflflll', 'lead', 'cry', 'dark', 'machine', 'note', 'waitwalan', 'fifife', 'star', 'box', 'noun', 'field', 'rest', 'correct', 'able', 'pound', 'done', 'beauty', 'drive', 'stood', 'contain', 'front', 'teach', 'week', 'final', 'gave', 'green', 'oh', 'quick', 'develop', 'ocean', 'warm', 'free', 'minute', 'strong', 'specispecisd', 'behind', 'cccccctail', 'produce', 'fact', 'street', 'inch', 'multiply', 'nothing', 'course', 'stay', 'wheel', 'full', 'force', 'blue', 'object', 'decide', 'surface', 'deep', 'moon', 'island', 'foot', 'system', 'busy', 'test', 'record', 'boat', 'common', 'gold', 'possible', 'plane', 'steasteay', 'wonder', 'laugh', 'thousand', 'ago', 'ran', 'check', 'game', 'shape', 'equate', 'hot', 'miss', 'brought', 'heat', 'snow', 'tire', 'bring', 'yes', 'distant', 'fififeast', 'paint', 'language', 'among', 'grand', 'ball', 'yet', 'yet', '', '', 'gop', 'heart', 'am', 'present', 'heaheadance', 'engine', 'position', 'arm', 'wide', 'sail', 'material', 'size', 'vary', 'settle', 'speak', 'weight', 'general', 'ice', 'matter', 'circle', 'pair', 'include', 'divide', 'syllable', 'felt', 'perhaps', 'pick', 'sudden', 'count', 'square', 'reason', 'length', 'represent', 'art', 'subject', 'region', 'energyenerg', 'probable', 'bed', 'brother', 'egg', 'ride', 'cell', 'believe', 'fraction', 'forest', 'sit', 'race', 'window', 'store', 'summer', 'train', 'sleep', 'prove', 'lone', 'lelelxercise', 'wall', 'catch', 'mount', 'wish', 'skyskyskd', 'joy', 'winter', 'sat', 'written', 'wild', 'instrument', 'kept', 'glass', 'grass', 'cow', 'job', 'edge', 'sign', 'visit', 'ppppppoft', 'fun', 'bright', 'gggggeather', 'month', 'million', 'bear', 'finish', 'happy', 'hope', 'flower', 'clothe', 'strange', 'gonegonmpgoney', 'eight', 'village', 'meet', 'root', 'buy', 'raise', 'solve', 'metal', 'whether', 'push', 'seven', 'paragraph', 'third', 'shall', 'held', 'hair', 'describe', 'cook', 'floor', 'either', 'result', 'burn', 'hill', 'safe', 'cat', 'century', 'consider', 'type', 'law', 'bit', 'coast', 'copy', 'phrase', 'silent', 'tall', 'sand', 'ssss', 'roll', 'temperature', 'ffffff', 'industry', 'value', 'fight', 'lie', 'beat', 'excite', 'naturalnaturalense', 'eee', 'else', 'ququq', 'bbbbb', 'case', 'middle', 'kill', 'son', 'lake', 'moment', 'scale', 'loud', 'spring', 'observe', 'child', 'straight', 'consonant', 'nation', 'dictionary', 'milk', 'speed', 'method', 'organ', 'pay', 'age', 'section', 'dress', 'cloud', 'surprsue', 'quiet', 'stone', 'tiny', 'climb', 'cool', 'design', 'ppppplot', 'experiment', 'bottom', 'key', 'iron', 'single', 'stick', 'flat', 'twenty', 'skin', 'smile', 'crease', 'hole', 'trade', 'melody', 'trip', 'office', 'receive', 'row', 'row', 'ive', 'act', 'symbol', 'die', 'least', 'trouble', 'shout', 'except', 'wrote', 'seed', 'tone', 'join', 'joigest', 'clean', 'break', 'lalalalrd', 'rise', 'badbadba', 'oil', 'blood', 'touch', 'grew', 'cent', 'mix', 'team', 'wire', 'cost', 'lost', 'brown', 'wear', 'garden', 'equal', 'sent', 'choose', 'fell', 'fit', 'flow', 'fair', 'bank', 'collect', 'save', 'control', 'decimal', 'gentle', 'woman', 'captain', 'practice', 'separatsepaffiseparatsepr', 'please', 'protect', 'noon', 'whose', 'locate', 'ring', 'character', 'insect', 'caught', 'period', 'indicate', 'radio', 'spoke', 'atom', 'human', 'history', 'effect', 'electric', 'expect', 'crop', 'modern', 'element', 'hit', 'student', 'corner', 'corner', 'upply', 'bone', 'rail', 'imagine', 'proproe', 'agree', 'thus', 'capital', "won't", 'chair', 'danger', 'fruit', 'rich', 'thick', 'thickerthickess', 'opeopee', 'guessguecessary', 'sharp', 'wing', 'create', 'neighbor', 'wash', 'bat', 'rather', 'crowd', 'corn', 'compare', 'poem', 'string', 'bell', 'depend', 'meat', 'rub', 'tube', 'famous', 'dollar', 'stream', 'fear', 'sight', 'thin', 'triangle', 'planet', 'hurry', 'chief', 'colony', 'clock', 'mine', 'tie', 'enter', 'major', 'fresh', 'search', 'send', 'yellow', 'gun', 'alloalloint', 'deaddeaddeaesert', 'suit', 'curcurt', 'lift', 'rose', 'continue', 'block', 'chart', 'hat', 'sell', 'succesu', 'company', 'subtrsubtevent', 'particular', 'deal', 'swim', 'term', 'opposite', 'wife', 'shoe', 'shoulder', 'spread', 'arrange', 'camp', 'invent', 'cotton', 'born', 'determine', 'quququnine', 'truck', 'noise', 'level', 'chance', 'chance', 'shop', 'stretch', 'throw', 'shine', 'property', 'column', 'molecule', 'selsel', 'wrong', 'gray', 'repeat', 'require', 'broad', 'prepare', 'salt', 'nose', 'plural', 'anger', 'claim', 'continent', 'oxygen', 'sugar', 'death', 'deatty', 'skill', 'women', 'season', 'solution', 'magnet', 'silver', 'thank', 'branch', 'match', 'suffix', 'especially', 'fig', 'afraid', 'huge', 'sister', 'steel', 'discuss', 'forward', 'similar', 'guide', 'experience', 'score', 'apple', 'bought', 'ledledled', 'coat', 'mass', 'card', 'babababpebabipbababdream', 'evening', 'condition', 'feed', 'tool', 'total', 'basic', 'smell', 'smell', '', 'nor', 'double', 'seat', 'arrive', 'master', 'track', 'parent', 'shore', 'division', 'sheet', 'substance', 'favor', 'connect', 'post', 'spend', 'chord', 'fat', 'glad', 'original', 'share', 'stationstad', 'bread', 'charge', 'proper', 'bar', 'offer', 'segmentsegave', 'duck', 'instant', 'market', 'degree', 'populate', 'chick', 'dear', 'enemy', 'reply', 'drink', 'occur', 'sssssrt', 'speech', 'nature', 'range', 'steam', 'motion', 'path', 'liquid', 'log', 'meant', 'quotient', 'teetteetteetteck']

        # Let's create a lot of objects!
        LOTS = 10000
        experiments = []
        xsa = []
        for x in range(1, LOTS):
            ex = Experiment()
            ex.accession_code = "".join(random.choice(words)
                                        for i in range(3)) + str(random.randint(0, 1000))[:64]
            ex.title = " ".join(random.choice(words) for i in range(10))
            ex.description = " ".join(random.choice(words) for i in range(100))
            ex.technology = random.choice(["RNA-SEQ", "MICROARRAY"])
            ex.submitter_institution = random.choice(["Funkytown", "Monkeytown"])
            experiments.append(ex)

        homo_sapiens = Organism.get_object_for_name("HOMO_SAPIENS")

        Experiment.objects.bulk_create(experiments)
        for exz in experiments:
            xs = ExperimentSampleAssociation()
            xs.experiment = exz
            xs.sample = sample
            xsa.append(xs)
        ExperimentSampleAssociation.objects.bulk_create(xsa)
        experiments = []
        xsa = []

        ex = Experiment()
        ex.accession_code = "FINDME_TEMPURA"
        ex.title = "THISWILLBEINASEARCHRESULT"
        ex.description = "SOWILLTHIS"
        ex.technology = "MICROARRAY"
        ex.submitter_institution = "Funkytown"
        experiments.append(ex)

        ex2 = Experiment()
        ex2.accession_code = "FINDME2"
        ex2.title = "THISWILLBEINASEARCHRESULT"
        ex2.description = "SOWILLTHIS"
        ex2.technology = "RNA-SEQ"
        ex2.submitter_institution = "Funkytown"
        ex2.has_publication = True
        experiments.append(ex2)

        ex3 = Experiment()
        ex3.accession_code = "FINDME3"
        ex3.title = "THISWILLBEINASEARCHRESULT"
        ex3.description = "SOWILLTHIS"
        ex3.technology = "FAKE-TECH"
        ex3.submitter_institution = "Utopia"
        experiments.append(ex3)

        # Use an E-GEOD-XXXX accession so we can test the E-GEOD -> GSE alternate accession.
        ex4 = Experiment()
        ex4.accession_code = "E-GEOD-1234"
        ex4.title = "IGNORED"
        ex4.description = "IGNORED"
        ex4.technology = "MICROARRAY"
        ex4.submitter_institution = "IGNORED"
        # Bulk create won't call the save method.
        ex4.save()

        # Use an GSEXXX accession so we can test the GSE -> E-GEOD alternate accession.
        ex5 = Experiment()
        ex5.accession_code = "GSE5678"
        ex5.title = "IGNORED"
        ex5.description = "IGNORED"
        ex5.technology = "RNA-SEQ"
        ex5.submitter_institution = "IGNORED"
        ex5.has_publication = True
        # Bulk create won't call the save method.
        ex5.save()

        sample1 = Sample()
        sample1.title = "1123"
        sample1.accession_code = "1123"
        sample1.platform_name = "AFFY"
        sample1.is_processed = True
        sample1.technology = "RNA-SEQ"
        sample1.save()

        sample2 = Sample()
        sample2.title = "3345"
        sample2.accession_code = "3345"
        sample2.platform_name = "ILLUMINA"
        sample2.organism = homo_sapiens
        sample2.is_processed = True
        sample1.technology = "MICROARRAY"
        sample2.save()

        Experiment.objects.bulk_create(experiments)
        experiment_sample_association = ExperimentSampleAssociation()
        experiment_sample_association.sample = sample
        experiment_sample_association.experiment = ex2
        experiment_sample_association.save()
        experiment_sample_association = ExperimentSampleAssociation()
        experiment_sample_association.sample = sample
        experiment_sample_association.experiment = ex3
        experiment_sample_association.save()
        experiment_sample_association = ExperimentSampleAssociation()
        experiment_sample_association.sample = sample
        experiment_sample_association.experiment = ex4
        experiment_sample_association.save()
        experiment_sample_association = ExperimentSampleAssociation()
        experiment_sample_association.sample = sample
        experiment_sample_association.experiment = ex5
        experiment_sample_association.save()

        xa = ExperimentAnnotation()
        xa.data = {'name': 'Clark Kent'}
        xa.experiment = ex
        xa.save()

        xoa = ExperimentOrganismAssociation()
        xoa.experiment=ex
        xoa.organism=homo_sapiens
        xoa.save()

        xoa = ExperimentOrganismAssociation()
        xoa.experiment=ex2
        xoa.organism=homo_sapiens
        xoa.save()

        xoa = ExperimentOrganismAssociation()
        xoa.experiment=ex3
        xoa.organism=Organism.objects.create(name="Extra-Terrestrial-1982", taxonomy_id=9999)
        xoa.save()

        experiment_sample_association = ExperimentSampleAssociation()
        experiment_sample_association.sample = sample1
        experiment_sample_association.experiment = ex
        experiment_sample_association.save()

        experiment_sample_association = ExperimentSampleAssociation()
        experiment_sample_association.sample = sample2
        experiment_sample_association.experiment = ex
        experiment_sample_association.save()

        # Test all
        response = self.client.get(reverse('search'))
        self.assertEqual(response.json()['count'], LOTS + 5)

        # Test search
        response = self.client.get(reverse('search'), {'search': 'THISWILLBEINASEARCHRESULT'})
        self.assertEqual(response.json()['count'], 3)

        # Test filter
        response = self.client.get(reverse('search'), {'search': 'FINDME', 'has_publication': True})
        for result in response.json()['results']:
            self.assertTrue(result['has_publication'])
        response = self.client.get(reverse('search'), {'search': 'FINDME', 'has_publication': False})
        for result in response.json()['results']:
            self.assertFalse(result['has_publication'])

        # Test search and filter
        response = self.client.get(reverse('search'),
                                   {'search': 'THISWILLBEINASEARCHRESULT',
                                    'technology': 'MICROARRAY'})
        self.assertEqual(response.json()['count'], 1)
        self.assertEqual(response.json()['results'][0]['accession_code'], 'FINDME_TEMPURA')
        # filter contain values for the top search query
        self.assertEqual(response.json()['filters']['technology'], {'FAKE-TECH': 1, 'MICROARRAY': 2, 'RNA-SEQ': 1})
        self.assertEqual(response.json()['filters']['publication'], {'has_publication': 1})
        self.assertEqual(response.json()['filters']['organism'], {'Extra-Terrestrial-1982': 1, 'HOMO_SAPIENS': 3})

        # Test search and interactive filtering
        response = self.client.get(reverse('search'),
                                   {'search': 'THISWILLBEINASEARCHRESULT',
                                    'technology': 'MICROARRAY',
                                    'filter_order': 'technology'})
        self.assertEqual(response.json()['count'], 1)
        self.assertEqual(response.json()['results'][0]['accession_code'], 'FINDME_TEMPURA')
        self.assertEqual(response.json()['filters']['technology'], {'FAKE-TECH': 1, 'MICROARRAY': 2, 'RNA-SEQ': 1})
        self.assertEqual(response.json()['filters']['publication'], {}) # FINDME_TEMPURA is the only result and doesn't have any publication
        # organism filter number should reflect the single result: two samples HOMO_SAPIENS
        self.assertEqual(response.json()['filters']['organism'], {'HOMO_SAPIENS': 2})

        response = self.client.get(reverse('search'),
                                   {'search': 'THISWILLBEINASEARCHRESULT',
                                    'organisms__name': 'Extra-Terrestrial-1982'})
        self.assertEqual(response.json()['count'], 1)

        response = self.client.get(reverse('search'), {'search': 'Clark Kent'})
        self.assertEqual(response.json()['count'], 1)
        self.assertEqual(response.json()['results'][0]['accession_code'], "FINDME_TEMPURA")

        # Test multiple filters
        # This has to be done manually due to dicts requring distinct keys
        response = self.client.get(reverse('search') + "?search=THISWILLBEINASEARCHRESULT&technology=MICROARRAY&technology=FAKE-TECH")
        self.assertEqual(response.json()['count'], 2)

        # Test ordering
        response = self.client.get(reverse('search') + "?search=SEARCH&ordering=id")
        response2 = self.client.get(reverse('search') + "?search=SEARCH&ordering=-id")
        self.assertNotEqual(response.json()['results'][0]['id'], response2.json()['results'][0]['id'])
        self.assertTrue(response2.json()['results'][0]['id'] > response.json()['results'][0]['id'])

        # Test Searching on Alternate Accession Codes
        response = self.client.get(reverse('search'), {'search': 'GSE1234'})
        self.assertEqual(response.json()['count'], 1)
        self.assertEqual(response.json()['results'][0]['accession_code'], "E-GEOD-1234")

        response = self.client.get(reverse('search'), {'search': 'E-GEOD-5678'})
        self.assertEqual(response.json()['count'], 1)
        self.assertEqual(response.json()['results'][0]['accession_code'], "GSE5678")

        cr = ComputationalResult()
        cr.save()

        qni = ComputedFile()
        qni.is_qn_target = True
        qni.s3_bucket = "fake_qni_bucket"
        qni.s3_key = "zazaza_homo_sapiens_1234.tsv"
        qni.filename = "homo_sapiens_1234.tsv"
        qni.is_public = True
        qni.size_in_bytes = 56789
        qni.sha1 = "c0a88d0bb020dadee3b707e647f7290368c235ba"
        qni.result = cr
        qni.save()

        qni = ComputedFile()
        qni.is_qn_target = False
        qni.s3_bucket = "X"
        qni.s3_key = "X.tsv"
        qni.filename = "XXXXXXXXXXXXXXX.tsv"
        qni.is_public = True
        qni.size_in_bytes = 1
        qni.sha1 = "123"
        qni.result = cr
        qni.save()

        response = self.client.get(reverse('qn-targets'))
        self.assertEqual(len(response.json()), 1)
        self.assertEqual(response.json()[0]['s3_url'], 'https://s3.amazonaws.com/fake_qni_bucket/zazaza_homo_sapiens_1234.tsv')

    @patch('data_refinery_common.message_queue.send_job')
    def test_create_update_dataset(self, mock_send_job):

        # Get a token first
        response = self.client.get(reverse('token'),
                                    content_type="application/json")
        token = response.json()
        token['is_activated'] = True
        token_id = token['id']
        response = self.client.post(reverse('token'),
                                    json.dumps(token),
                                    content_type="application/json")

        activated_token = response.json()
        self.assertEqual(activated_token['id'], token_id)
        self.assertEqual(activated_token['is_activated'], True)

        # Good
        jdata = json.dumps({'data': {"A": ["B"]}})
        response = self.client.post(reverse('create_dataset'),
                                    jdata,
                                    content_type="application/json")

        self.assertEqual(response.status_code, 201)
        self.assertEqual(response.json()['data'], json.loads(jdata)['data'])
        good_id = response.json()['id']

        response = self.client.get(reverse('dataset', kwargs={'id': good_id}))
        self.assertEqual(response.json()['id'], good_id)
        self.assertEqual(response.json()['data'], json.loads(jdata)['data'])
        self.assertEqual(response.json()['data']["A"], ["B"])

        # Bad (Duplicates)
        jdata = json.dumps({'data': {"A": ["B", "B", "B"]}})
        response = self.client.post(reverse('create_dataset'),
                                    jdata,
                                    content_type="application/json")

        self.assertEqual(response.status_code, 400)

        # Update, just an experiment accession
        jdata = json.dumps({'data': {"GSE123": ["ALL"]}})
        response = self.client.put(reverse('dataset', kwargs={'id': good_id}),
                                   jdata,
                                   content_type="application/json")
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.json()['id'], good_id)
        # We are not entirely RESTful here, that's okay.
        self.assertNotEqual(response.json()['data'], json.loads(jdata)['data'])
        self.assertEqual(response.json()['data']["GSE123"], ["789"])

        # Update
        jdata = json.dumps({'data': {"A": ["C"]}})
        response = self.client.put(reverse('dataset', kwargs={'id': good_id}),
                                   jdata,
                                   content_type="application/json")
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.json()['id'], good_id)
        self.assertEqual(response.json()['data'], json.loads(jdata)['data'])
        self.assertEqual(response.json()['data']["A"], ["C"])

        # Can't update if started
        dataset = Dataset.objects.get(id=good_id)
        dataset.is_processing = True
        dataset.save()
        jdata = json.dumps({'data': {"A": ["D"]}})
        response = self.client.put(reverse('dataset', kwargs={'id': good_id}),
                                   jdata,
                                   content_type="application/json")
        self.assertNotEqual(response.json()['data']["A"], ["D"])

        # Bad
        jdata = json.dumps({'data': 123})
        response = self.client.post(reverse('create_dataset'),
                                    jdata,
                                    content_type="application/json")
        self.assertEqual(response.status_code, 400)

        # This will actually kick off a job if we don't patch send_job or supply no_send_job
        dataset = Dataset.objects.get(id=good_id)
        dataset.is_processing = False
        dataset.save()

        # With bad token first
        jdata = json.dumps({'data': {"A": ["D"]}, 'start': True, 'no_send_job': True, 'token_id': "HEYO" } )
        response = self.client.put(reverse('dataset', kwargs={'id': good_id}), jdata, content_type="application/json")
        self.assertEqual(response.status_code, 500)

        jdata = json.dumps({'data': {"A": ["D"]}, 'start': True, 'no_send_job': True, 'token_id': token_id, 'email_address': 'trust@verify.com' } )
        response = self.client.put(reverse('dataset', kwargs={'id': good_id}), jdata, content_type="application/json")
        self.assertEqual(response.json()["is_processing"], True)

        ds = Dataset.objects.get(id=response.json()['id'])
        self.assertEqual(ds.email_address, 'trust@verify.com')

    def test_processed_samples_only(self):
        """ Don't return unprocessed samples """
        experiment = Experiment()
        experiment.accession_code = "GSX12345"
        experiment.is_public = True
        experiment.save()

        sample = Sample()
        sample.title = "I am unprocessed"
        sample.accession_code = "GSXUnprocessed"
        sample.is_processed = False
        sample.save()

        experiment_sample_association = ExperimentSampleAssociation()
        experiment_sample_association.sample = sample
        experiment_sample_association.experiment = experiment
        experiment_sample_association.save()

        response = self.client.get(reverse('search'), {'search': "GSX12345"})
        self.assertEqual(response.json()['count'], 0)

        sample2 = Sample()
        sample2.title = "I am processed"
        sample2.accession_code = "GSXProcessed"
        sample2.is_processed = True
        sample2.save()

        experiment_sample2_association = ExperimentSampleAssociation()
        experiment_sample2_association.sample = sample2
        experiment_sample2_association.experiment = experiment
        experiment_sample2_association.save()

        response = self.client.get(reverse('search'), {'search': "GSX12345"})
        self.assertEqual(response.json()['count'], 1)

        qs = Experiment.processed_public_objects
        self.assertEqual(len(experiment.processed_samples), 1)

        experiment.delete()
        sample.delete()
        sample2.delete()

    def test_dataset_stats(self):
        """ Test the dataset stats endpoint """

        homo_sapiens = Organism.get_object_for_name("HOMO_SAPIENS")
        time.sleep(30)
        gallus_gallus = Organism.get_object_for_name("GALLUS_GALLUS")
        time.sleep(30)
        equus_ferus = Organism.get_object_for_name("EQUUS_FERUS")
        time.sleep(30)

        ex = Experiment()
        ex.accession_code = "XYZ123"
        ex.title = "XYZ123"
        ex.description = "XYZ123"
        ex.technology = "MICROARRAY"
        ex.submitter_institution = "XYZ123"
        ex.save()

        ex2 = Experiment()
        ex2.accession_code = "ABC789"
        ex2.title = "ABC789"
        ex2.description = "ABC789"
        ex2.technology = "RNA-SEQ"
        ex2.submitter_institution = "Funkytown"
        ex2.save()

        sample1 = Sample()
        sample1.title = "1"
        sample1.accession_code = "1"
        sample1.platform_name = "AFFY"
        sample1.organism = homo_sapiens
        sample1.save()

        sample2 = Sample()
        sample2.title = "2"
        sample2.accession_code = "2"
        sample2.platform_name = "ILLUMINA"
        sample2.organism = gallus_gallus
        sample2.save()

        sample3 = Sample()
        sample3.title = "3"
        sample3.accession_code = "3"
        sample3.platform_name = "ILLUMINA"
        sample3.organism = gallus_gallus
        sample3.save()

        xoa = ExperimentOrganismAssociation()
        xoa.experiment=ex
        xoa.organism=homo_sapiens
        xoa.save()

        xoa = ExperimentOrganismAssociation()
        xoa.experiment=ex2
        xoa.organism=gallus_gallus
        xoa.save()

        xoa = ExperimentOrganismAssociation()
        xoa.experiment=ex2
        xoa.organism=equus_ferus
        xoa.save()

        experiment_sample_association = ExperimentSampleAssociation()
        experiment_sample_association.sample = sample1
        experiment_sample_association.experiment = ex
        experiment_sample_association.save()

        experiment_sample_association = ExperimentSampleAssociation()
        experiment_sample_association.sample = sample2
        experiment_sample_association.experiment = ex2
        experiment_sample_association.save()

        experiment_sample_association = ExperimentSampleAssociation()
        experiment_sample_association.sample = sample3
        experiment_sample_association.experiment = ex2
        experiment_sample_association.save()

        jdata = json.dumps({'data': {"XYZ123": ["1"], "ABC789": ["2"]}})
        response = self.client.post(reverse('create_dataset'),
                                    jdata,
                                    content_type="application/json")

        self.assertEqual(response.status_code, 201)
        self.assertEqual(response.json()['data'], json.loads(jdata)['data'])
        good_id = response.json()['id']

        response = self.client.get(reverse('dataset_stats', kwargs={'id': good_id}))
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        self.assertEqual(response.json()['GALLUS_GALLUS'], {'num_experiments': 1, 'num_samples': 1})
        self.assertEqual(response.json()['HOMO_SAPIENS'], {'num_experiments': 1, 'num_samples': 1})
        self.assertEqual(len(response.json().keys()), 2)

        # Check that we can fetch these sample details via samples API
        response = self.client.get(reverse('samples'), {'dataset_id': good_id})
        self.assertEqual(response.json()['count'], 2)

    @patch('raven.contrib.django.models.client')
    def test_sentry_middleware_ok(self, mock_client):
        # We don't even import raven if it's a good response.
        response = self.client.get(reverse('experiments'))
        self.assertEqual(response.status_code, status.HTTP_200_OK)
        mock_client.is_enabled.assert_not_called()

    @patch('raven.contrib.django.models.client')
    def test_sentry_middleware_404(self, mock_client):
        # We don't send anything to raven if it's not enabled
        mock_client.is_enabled.side_effect = lambda: False
        response = self.client.get(reverse('experiments_detail', kwargs={'pk': '1000'}))
        self.assertEqual(response.status_code, 404)
        mock_client.captureMessage.assert_not_called()

        # A 404 with raven enabled will send a message to sentry
        mock_client.is_enabled.side_effect = lambda: True
        response = self.client.get(reverse('experiments_detail', kwargs={'pk': '1000'}))
        self.assertEqual(response.status_code, 404)
        mock_client.captureMessage.assert_not_called()

        # A 404 with raven enabled will send a message to sentry
        mock_client.is_enabled.side_effect = lambda: True
        response = self.client.get(reverse('experiments_detail', kwargs={'pk': '1000'})[:-1] + "aasdas/")
        self.assertEqual(response.status_code, 404)
        mock_client.captureMessage.assert_not_called()

    @patch.object(ExperimentList, 'get')
    @patch('raven.contrib.django.models.client')
    def test_sentry_middleware_403(self, mock_client, mock_get_method):
        mock_get_method.side_effect = lambda _: HttpResponseForbidden()
        # A 403 with raven enabled will send a message to sentry
        mock_client.is_enabled.side_effect = lambda: True
        response = self.client.get(reverse('experiments'))
        self.assertEqual(response.status_code, 403)
        mock_client.captureMessage.assert_called()

    @patch.object(ExperimentList, 'get')
    @patch('raven.contrib.django.models.client')
    def test_sentry_middleware_500(self, mock_client, mock_get_method):
        def raise_error(_):
            raise KeyError()

        mock_get_method.side_effect = lambda _: HttpResponseServerError()
        # A 500 with raven enabled will send a message to sentry
        mock_client.is_enabled.side_effect = lambda: True
        response = self.client.get(reverse('experiments'))
        self.assertEqual(response.status_code, 500)
        mock_client.captureMessage.assert_called()

class ESTestCases(APITestCase):

    def test_es_endpoint(self):
        """ Test basic ES functionality 

        This is pretty tricky because ES doesn't know that we're creating
        test objects.
        """

        # First, we purge.
        Experiment.objects.all().delete()

        experiment = Experiment()
        experiment.accession_code = "GSE000-X"
        experiment.title = "NONONONO"
        experiment.description = "Boooooourns. Wasabi."
        experiment.technology = "RNA-SEQ"
        experiment.save()

        experiment = Experiment()
        experiment.accession_code = "GSE123-X"
        experiment.title = "Hey Ho Let's Go"
        experiment.description = "This is a very exciting test experiment. Faygo soda. Blah blah blah."
        experiment.technology = "MICROARRAY"
        experiment.save()
        self.experiment = experiment

        experiment_annotation = ExperimentAnnotation()
        experiment_annotation.data = {"hello": "world", "123": 456}
        experiment_annotation.experiment = experiment
        experiment_annotation.save()

        sample = Sample()
        sample.title = "123"
        sample.accession_code = "123"
        sample.save()

        sample = Sample()
        sample.title = "789"
        sample.accession_code = "789"
        sample.save()
        self.sample = sample

        sample_annotation = SampleAnnotation()
        sample_annotation.data = {"goodbye": "world", "789": 123}
        sample_annotation.sample = sample
        sample_annotation.save()

        original_file = OriginalFile()
        original_file.save()

        original_file_sample_association = OriginalFileSampleAssociation()
        original_file_sample_association.sample = sample
        original_file_sample_association.original_file = original_file
        original_file_sample_association.save()

        downloader_job = DownloaderJob()
        downloader_job.save()

        download_assoc = DownloaderJobOriginalFileAssociation()
        download_assoc.original_file = original_file
        download_assoc.downloader_job = downloader_job
        download_assoc.save()

        processor_job = ProcessorJob()
        processor_job.save()

        processor_assoc = ProcessorJobOriginalFileAssociation()
        processor_assoc.original_file = original_file
        processor_assoc.processor_job = processor_job
        processor_assoc.save()

        experiment_sample_association = ExperimentSampleAssociation()
        experiment_sample_association.sample = sample
        experiment_sample_association.experiment = experiment
        experiment_sample_association.save()

        result = ComputationalResult()
        result.save()

        sra = SampleResultAssociation()
        sra.sample = sample
        sra.result = result
        sra.save()

        result = ComputationalResult()
        result.save()

        sra = SampleResultAssociation()
        sra.sample = sample
        sra.result = result
        sra.save()

        # TODO: Use `reverse` when not using the DefaultRouter
        response = self.client.get('/es/')

        """ Test basic ES functionality """
        es_search_result = ExperimentDocument.search().filter("term", description="soda")
        es_search_result_qs = es_search_result.to_queryset()
        self.assertEqual(len(es_search_result_qs), 1)

        # Sanity
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.json()['count'], 2)

        # Basic Search
        response = self.client.get('/es/?search=soda')
        self.assertEqual(response.json()['count'], 1)

        # Positive filter result
        response = self.client.get('/es/?search=soda&technology=MICROARRAY')
        self.assertEqual(response.json()['count'], 1)

        # Negative filter result
        response = self.client.get('/es/?search=soda&technology=rna')
        self.assertEqual(response.json()['count'], 0)

class ProcessorTestCases(APITestCase):
    def setUp(self):
        salmon_quant_env = {
            'os_distribution': 'Ubuntu 16.04.4 LTS',
            'os_pkg': {
                'python3': '3.5.1-3',
                'python3-pip': '8.1.1-2ubuntu0.4'
            },
            'cmd_line': {
                'salmon --version': 'salmon 0.9.1'
            },
            'python': {
                'Django': '2.0.6',
                'data-refinery-common': '0.5.0'
            }
        }
        self.salmon_quant_proc = Processor.objects.create(
            name="Salmon Quant",
            version="0.45",
            docker_image="ccdl/salmon_img:v1.23",
            environment=salmon_quant_env
        )

        salmontools_env = {
            'os_distribution': 'Ubuntu 16.04.4 LTS',
            'os_pkg': {
                'python3': '3.5.1-3',
                'python3-pip': '8.1.1-2ubuntu0.4',
                'g++': '4:5.3.1-1ubuntu1',
                'cmake': '3.5.1-1ubuntu3'
            },
            'cmd_line': {
                'salmontools --version': 'Salmon Tools 0.1.0'
            },
            'python': {
                'Django': '2.0.6',
                'data-refinery-common': '0.5.0'
            }
        }
        Processor.objects.create(
            name="Salmontools",
            version="1.83",
            docker_image="ccdl/salmontools_img:v0.45",
            environment=salmontools_env
        )

    def test_endpoint(self):
        response = self.client.get(reverse('processors'))
        self.assertEqual(response.status_code, status.HTTP_200_OK)

        processors = response.json()
        self.assertEqual(processors[0]['name'], 'Salmon Quant')
        self.assertEqual(processors[0]['environment']['os_pkg']['python3'], '3.5.1-3')

        self.assertEqual(processors[1]['name'], 'Salmontools')
        self.assertEqual(processors[1]['environment']['cmd_line']['salmontools --version'],
                         'Salmon Tools 0.1.0')

    def test_processor_and_organism_in_sample(self):
        sample = Sample.objects.create(title="fake sample")
        organism = Organism.get_object_for_name("HOMO_SAPIENS")
        transcriptome_result = ComputationalResult.objects.create()
        organism_index = OrganismIndex.objects.create(organism=organism,
                                                      result=transcriptome_result,
                                                      index_type="TRANSCRIPTOME_LONG")
        result = ComputationalResult.objects.create(processor=self.salmon_quant_proc,
                                                    organism_index=organism_index)
        sra = SampleResultAssociation.objects.create(sample=sample, result=result)

        response = self.client.get(reverse('samples_detail',
                                           kwargs={'pk': sample.id}))
        self.assertEqual(response.status_code, status.HTTP_200_OK)

        processor = response.json()['results'][0]['processor']
        self.assertEqual(processor['name'], self.salmon_quant_proc.name)
        self.assertEqual(processor['environment']['os_pkg']['python3'],
                         self.salmon_quant_proc.environment['os_pkg']['python3'])

        organism_index = response.json()['results'][0]['organism_index']
        self.assertEqual(organism_index["result"], transcriptome_result.id)
        self.assertEqual(organism_index["index_type"], "TRANSCRIPTOME_LONG")
