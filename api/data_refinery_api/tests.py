from django.contrib.auth.models import User
from django.urls import reverse
from rest_framework import status
from rest_framework.test import APITestCase
from unittest.mock import patch

import json
import random

from data_refinery_common.models import (
    Experiment,
    ExperimentAnnotation,
    Sample,
    SampleAnnotation,
    ExperimentSampleAssociation,
    Organism,
    OriginalFile,
    OriginalFileSampleAssociation,
    DownloaderJob,
    DownloaderJobOriginalFileAssociation,
    ProcessorJob,
    ProcessorJobOriginalFileAssociation,
    ComputationalResult,
    Dataset
)
from data_refinery_api.serializers import (
    ExperimentSerializer,
    DetailedExperimentSerializer,
    SampleSerializer,
    DetailedSampleSerializer,
    OrganismSerializer,
    PlatformSerializer,
    InstitutionSerializer,

    # Jobs
    SurveyJobSerializer,
    DownloaderJobSerializer,
    ProcessorJobSerializer
)

class APITestCases(APITestCase):
    def setUp(self):
        # Saving this for if we have protected endpoints
        # self.superuser = User.objects.create_superuser('john', 'john@snow.com', 'johnpassword')
        # self.client.login(username='john', password='johnpassword')
        # self.user = User.objects.create(username="mike")

        experiment = Experiment()
        experiment.save()

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

        return

    def test_all_endpoints(self):
        response = self.client.get(reverse('experiments'))
        self.assertEqual(response.status_code, status.HTTP_200_OK)

        response = self.client.get(reverse('experiments'), kwargs={'page': 1})
        self.assertEqual(response.status_code, status.HTTP_200_OK)

        response = self.client.get(reverse('experiments_detail', kwargs={'pk': '1'}))
        self.assertEqual(response.status_code, status.HTTP_200_OK)

        response = self.client.get(reverse('samples'))
        self.assertEqual(response.status_code, status.HTTP_200_OK)

        response = self.client.get(reverse('samples'), {'ids': str(self.sample.id) + ',1000'})
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


    def test_search_and_filter(self):

        # Our Docker image doesn't have the standard dict. >=[
        words = ['the', 'of', 'to', 'and', 'a', 'in', 'is', 'it', 'you', 'that', 'he', 'was', 'for', 'on', 'are', 'with', 'as', 'I', 'his', 'they', 'be', 'at', 'one', 'have', 'this', 'from', 'or', 'had', 'by', 'hot', 'word', 'but', 'what', 'some', 'we', 'can', 'out', 'other', 'were', 'all', 'there', 'when', 'up', 'use', 'your', 'how', 'said', 'an', 'each', 'she', 'which', 'do', 'their', 'time', 'if', 'will', 'way', 'about', 'many', 'then', 'them', 'write', 'would', 'like', 'so', 'these', 'her', 'long', 'make', 'thing', 'see', 'him', 'two', 'has', 'look', 'more', 'day', 'could', 'go', 'come', 'did', 'number', 'sound', 'no', 'most', 'people', 'my', 'over', 'know', 'water', 'than', 'call', 'first', 'who', 'may', 'down', 'side', 'been', 'now', 'find', 'any', 'new', 'work', 'part', 'take', 'get', 'place', 'made', 'live', 'where', 'after', 'back', 'little', 'only', 'round', 'man', 'year', 'came', 'show', 'every', 'good', 'me', 'give', 'our', 'under', 'name', 'very', 'through', 'just', 'form', 'sentence', 'great', 'think', 'say', 'help', 'low', 'line', 'differ', 'turn', 'cause', 'much', 'mean', 'before', 'move', 'right', 'boy', 'old', 'too', 'same', 'tell', 'does', 'set', 'three', 'want', 'air', 'well', 'also', 'play', 'small', 'end', 'put', 'home', 'read', 'hand', 'port', 'large', 'spell', 'add', 'even', 'land', 'here', 'must', 'big', 'high', 'such', 'follow', 'act', 'why', 'ask', 'men', 'change', 'went', 'light', 'kind', 'off', 'need', 'house', 'picture', 'try', 'us', 'again', 'animal', 'point', 'mother', 'world', 'near', 'build', 'self', 'earth', 'father', 'head', 'stand', 'own', 'page', 'should', 'country', 'found', 'answer', 'school', 'grow', 'study', 'still', 'learn', 'plant', 'cover', 'food', 'sun', 'four', 'between', 'state', 'keep', 'eye', 'never', 'last', 'let', 'thought', 'citcitree', 'cross', 'farm', 'hhhh', 'start', 'might', 'stosy', 'saw', 'sar', 'sea', 'draw', 'left', 'late', 'run', "don't", 'while', 'press', 'close', 'night', 'real', 'life', 'few', 'north', 'open', 'seemseegether', 'next', 'white', 'chilchin', 'chiln', 'got', 'walk', 'exampexamase', 'paper', 'group', 'always', 'music', 'those', 'botbotark', 'often', 'letter', 'until', 'mile', 'river', 'car', 'feet', 'care', 'second', 'book', 'carry', 'took', 'science', 'eat', 'room', 'friefr', 'bbban', 'idea', 'fish', 'mountain', 'stop', 'once', 'base', 'hear', 'horse', 'cut', 'sure', 'watch', 'color', 'face', 'wood', 'main', 'enough', 'plain', 'girl', 'usual', 'young', 'ready', 'above', 'ever', 'red', 'list', 'though', 'feel', 'talk', 'bird', 'soon', 'body', 'dog', 'family', 'direct', 'pose', 'leave', 'song', 'measure', 'door', 'product', 'black', 'short', 'numeral', 'class', 'wind', 'question', 'happen', 'complete', 'ship', 'area', 'half', 'rock', 'order', 'fire', 'south', 'problem', 'piece', 'told', 'knew', 'pass', 'since', 'top', 'whole', 'king', 'space', 'heard', 'best', 'hour', 'better', 'true', 'during', 'hundred', 'five', 'rrrember', 'step', 'early', 'hold', 'west', 'groundgroterest', 'reach', 'fast', 'verb', 'sing', 'llsten', 'six', 'table', 'travel', 'less', 'morning', 'ten', 'simple', 'several', 'vowel', 'toward', 'war', 'lay', 'against', 'pattern', 'slow', 'center', 'love', 'person', 'money', 'serve', 'appear', 'road', 'map', 'rain', 'rule', 'govern', 'pull', 'cold', 'notice', 'voice', 'unit', 'powepotown', 'fine', 'certain', 'flflflll', 'lead', 'cry', 'dark', 'machine', 'note', 'waitwalan', 'fifife', 'star', 'box', 'noun', 'field', 'rest', 'correct', 'able', 'pound', 'done', 'beauty', 'drive', 'stood', 'contain', 'front', 'teach', 'week', 'final', 'gave', 'green', 'oh', 'quick', 'develop', 'ocean', 'warm', 'free', 'minute', 'strong', 'specispecisd', 'behind', 'cccccctail', 'produce', 'fact', 'street', 'inch', 'multiply', 'nothing', 'course', 'stay', 'wheel', 'full', 'force', 'blue', 'object', 'decide', 'surface', 'deep', 'moon', 'island', 'foot', 'system', 'busy', 'test', 'record', 'boat', 'common', 'gold', 'possible', 'plane', 'steasteay', 'wonder', 'laugh', 'thousand', 'ago', 'ran', 'check', 'game', 'shape', 'equate', 'hot', 'miss', 'brought', 'heat', 'snow', 'tire', 'bring', 'yes', 'distant', 'fififeast', 'paint', 'language', 'among', 'grand', 'ball', 'yet', 'yet', '', '', 'gop', 'heart', 'am', 'present', 'heaheadance', 'engine', 'position', 'arm', 'wide', 'sail', 'material', 'size', 'vary', 'settle', 'speak', 'weight', 'general', 'ice', 'matter', 'circle', 'pair', 'include', 'divide', 'syllable', 'felt', 'perhaps', 'pick', 'sudden', 'count', 'square', 'reason', 'length', 'represent', 'art', 'subject', 'region', 'energyenerg', 'probable', 'bed', 'brother', 'egg', 'ride', 'cell', 'believe', 'fraction', 'forest', 'sit', 'race', 'window', 'store', 'summer', 'train', 'sleep', 'prove', 'lone', 'lelelxercise', 'wall', 'catch', 'mount', 'wish', 'skyskyskd', 'joy', 'winter', 'sat', 'written', 'wild', 'instrument', 'kept', 'glass', 'grass', 'cow', 'job', 'edge', 'sign', 'visit', 'ppppppoft', 'fun', 'bright', 'gggggeather', 'month', 'million', 'bear', 'finish', 'happy', 'hope', 'flower', 'clothe', 'strange', 'gonegonmpgoney', 'eight', 'village', 'meet', 'root', 'buy', 'raise', 'solve', 'metal', 'whether', 'push', 'seven', 'paragraph', 'third', 'shall', 'held', 'hair', 'describe', 'cook', 'floor', 'either', 'result', 'burn', 'hill', 'safe', 'cat', 'century', 'consider', 'type', 'law', 'bit', 'coast', 'copy', 'phrase', 'silent', 'tall', 'sand', 'ssss', 'roll', 'temperature', 'ffffff', 'industry', 'value', 'fight', 'lie', 'beat', 'excite', 'naturalnaturalense', 'eee', 'else', 'ququq', 'bbbbb', 'case', 'middle', 'kill', 'son', 'lake', 'moment', 'scale', 'loud', 'spring', 'observe', 'child', 'straight', 'consonant', 'nation', 'dictionary', 'milk', 'speed', 'method', 'organ', 'pay', 'age', 'section', 'dress', 'cloud', 'surprsue', 'quiet', 'stone', 'tiny', 'climb', 'cool', 'design', 'ppppplot', 'experiment', 'bottom', 'key', 'iron', 'single', 'stick', 'flat', 'twenty', 'skin', 'smile', 'crease', 'hole', 'trade', 'melody', 'trip', 'office', 'receive', 'row', 'row', 'ive', 'act', 'symbol', 'die', 'least', 'trouble', 'shout', 'except', 'wrote', 'seed', 'tone', 'join', 'joigest', 'clean', 'break', 'lalalalrd', 'rise', 'badbadba', 'oil', 'blood', 'touch', 'grew', 'cent', 'mix', 'team', 'wire', 'cost', 'lost', 'brown', 'wear', 'garden', 'equal', 'sent', 'choose', 'fell', 'fit', 'flow', 'fair', 'bank', 'collect', 'save', 'control', 'decimal', 'gentle', 'woman', 'captain', 'practice', 'separatsepaffiseparatsepr', 'please', 'protect', 'noon', 'whose', 'locate', 'ring', 'character', 'insect', 'caught', 'period', 'indicate', 'radio', 'spoke', 'atom', 'human', 'history', 'effect', 'electric', 'expect', 'crop', 'modern', 'element', 'hit', 'student', 'corner', 'corner', 'upply', 'bone', 'rail', 'imagine', 'proproe', 'agree', 'thus', 'capital', "won't", 'chair', 'danger', 'fruit', 'rich', 'thick', 'thickerthickess', 'opeopee', 'guessguecessary', 'sharp', 'wing', 'create', 'neighbor', 'wash', 'bat', 'rather', 'crowd', 'corn', 'compare', 'poem', 'string', 'bell', 'depend', 'meat', 'rub', 'tube', 'famous', 'dollar', 'stream', 'fear', 'sight', 'thin', 'triangle', 'planet', 'hurry', 'chief', 'colony', 'clock', 'mine', 'tie', 'enter', 'major', 'fresh', 'search', 'send', 'yellow', 'gun', 'alloalloint', 'deaddeaddeaesert', 'suit', 'curcurt', 'lift', 'rose', 'continue', 'block', 'chart', 'hat', 'sell', 'succesu', 'company', 'subtrsubtevent', 'particular', 'deal', 'swim', 'term', 'opposite', 'wife', 'shoe', 'shoulder', 'spread', 'arrange', 'camp', 'invent', 'cotton', 'born', 'determine', 'quququnine', 'truck', 'noise', 'level', 'chance', 'chance', 'shop', 'stretch', 'throw', 'shine', 'property', 'column', 'molecule', 'selsel', 'wrong', 'gray', 'repeat', 'require', 'broad', 'prepare', 'salt', 'nose', 'plural', 'anger', 'claim', 'continent', 'oxygen', 'sugar', 'death', 'deatty', 'skill', 'women', 'season', 'solution', 'magnet', 'silver', 'thank', 'branch', 'match', 'suffix', 'especially', 'fig', 'afraid', 'huge', 'sister', 'steel', 'discuss', 'forward', 'similar', 'guide', 'experience', 'score', 'apple', 'bought', 'ledledled', 'coat', 'mass', 'card', 'babababpebabipbababdream', 'evening', 'condition', 'feed', 'tool', 'total', 'basic', 'smell', 'smell', '', 'nor', 'double', 'seat', 'arrive', 'master', 'track', 'parent', 'shore', 'division', 'sheet', 'substance', 'favor', 'connect', 'post', 'spend', 'chord', 'fat', 'glad', 'original', 'share', 'stationstad', 'bread', 'charge', 'proper', 'bar', 'offer', 'segmentsegave', 'duck', 'instant', 'market', 'degree', 'populate', 'chick', 'dear', 'enemy', 'reply', 'drink', 'occur', 'sssssrt', 'speech', 'nature', 'range', 'steam', 'motion', 'path', 'liquid', 'log', 'meant', 'quotient', 'teetteetteetteck']

        # Let's create a lot of objects!
        LOTS = 10000
        experiments = []
        for x in range(1, LOTS):
            ex = Experiment()
            ex.accession_code = "".join(random.choice(words) for i in range(3)) + str(random.randint(0, 1000))[:64]
            ex.title =  " ".join(random.choice(words) for i in range(10))
            ex.description = " ".join(random.choice(words) for i in range(100))
            ex.technology = random.choice(["RNA-SEQ", "MICROARRAY"])
            ex.submitter_institution = random.choice(["Funkytown", "Monkeytown"])
            experiments.append(ex)

        ex = Experiment()
        ex.accession_code = "FINDME"
        ex.title =  "THISWILLBEINASEARCHRESULT"
        ex.description =  "SOWILLTHIS"
        ex.technology ="MICROARRAY"
        ex.submitter_institution = "Funkytown"
        experiments.append(ex)

        ex = Experiment()
        ex.accession_code = "FINDME2"
        ex.title =  "THISWILLBEINASEARCHRESULT"
        ex.description =  "SOWILLTHIS"
        ex.technology = "RNA-SEQ"
        ex.submitter_institution = "Funkytown"
        experiments.append(ex)

        Experiment.objects.bulk_create(experiments)

        # Test all
        response = self.client.get(reverse('search'))
        self.assertEqual(response.json()['count'], LOTS + 2)

        # Test search
        response = self.client.get(reverse('search'), {'search': 'THISWILLBEINASEARCHRESULT'})
        self.assertEqual(response.json()['count'], 2)

        # Test search and filter
        response = self.client.get(reverse('search'), {'search': 'THISWILLBEINASEARCHRESULT', 'technology': 'MICROARRAY'})
        self.assertEqual(response.json()['count'], 1)
        self.assertEqual(response.json()['results'][0]['accession_code'], 'FINDME')

    @patch('data_refinery_common.message_queue.send_job')
    def test_create_update_dataset(self, mock_send_job):

        # Good
        jdata = json.dumps({'data': {"A": ["B"]}})
        response = self.client.post(reverse('create_dataset'), jdata, content_type="application/json")

        self.assertEqual(response.status_code, 201)
        self.assertEqual(response.json()['data'], json.loads(jdata)['data'])
        good_id = response.json()['id']

        response = self.client.get(reverse('dataset', kwargs={'id': good_id}))
        self.assertEqual(response.json()['id'], good_id)
        self.assertEqual(response.json()['data'], json.loads(jdata)['data'])
        self.assertEqual(response.json()['data']["A"], ["B"])

        # Update
        jdata = json.dumps({'data': {"A": ["C"]}})
        response = self.client.put(reverse('dataset', kwargs={'id': good_id}), jdata, content_type="application/json")
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.json()['id'], good_id)
        self.assertEqual(response.json()['data'], json.loads(jdata)['data'])
        self.assertEqual(response.json()['data']["A"], ["C"])

        # Can't update if started
        dataset = Dataset.objects.get(id=good_id)
        dataset.is_processing = True
        dataset.save()
        jdata = json.dumps({'data': {"A": ["D"]}})
        response = self.client.put(reverse('dataset', kwargs={'id': good_id}), jdata, content_type="application/json")
        self.assertNotEqual(response.json()['data']["A"], ["D"])

        # Bad
        jdata = json.dumps({'data': 123})
        response = self.client.post(reverse('create_dataset'), jdata, content_type="application/json")
        self.assertEqual(response.status_code, 400)

        # This will actually kick off a job if we don't patch send_job
        dataset = Dataset.objects.get(id=good_id)
        dataset.is_processing = False
        dataset.save()
        jdata = json.dumps({'data': {"A": ["D"]}, 'start': True, 'no_send_job': True} )
        response = self.client.put(reverse('dataset', kwargs={'id': good_id}), jdata, content_type="application/json")
        self.assertEqual(response.json()["is_processing"], True)

