import React, { Component } from 'react';
import ChartistGraph from 'react-chartist';
import './App.css';

class App extends Component {
  componentWillMount() {
    this.fetchQueues();
  }

  fetchQueues = async () => {
    try {
      const res = await fetch('/v1/evaluations');
      console.log(res);
    } catch (e) {
      console.log(e);
    }
  };

  render() {
    const dataBar = {
      labels: ['Queue'],
      series: [[800000], [200000], [100000]]
    };

    const optionsBar = {
      stackBars: true
    };
    return (
      <div className="App">
        <header className="App-header">
          <h1 className="App-title">Executive Dashboard</h1>
        </header>
        <ChartistGraph data={dataBar} type="Bar" options={optionsBar} />
      </div>
    );
  }
}

export default App;
