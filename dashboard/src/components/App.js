import React, { Component } from 'react';
import './App.scss';
import './custom-chartist.scss';
import { ResponsiveContainer, PieChart, Pie, Cell, Tooltip } from 'recharts';

class App extends Component {
  constructor(props) {
    super(props);
    this.state = {
      totalQueueLength: []
    };
  }

  async componentDidMount() {
    this.setState({
      totalQueueLength: await this.fetchQueues()
    });
  }

  fetchQueues = async () => {
    try {
      const res = await (await fetch('../stats.json')).json();
      const totalQueueLengthData = Object.keys(res).map((job, i) => {
        return {
          name: job,
          value: res[job].total
        };
      });
      return totalQueueLengthData;
    } catch (e) {
      console.log(e);
    }
  };

  render() {
    const COLORS = ['#0088FE', '#00C49F', '#FFBB28', '#FF8042'];

    return (
      <div className="App">
        <header className="App-header">
          <h1 className="App-title">Executive Dashboard</h1>
        </header>
        <h2>Total Length of Queues by Type</h2>
        <div
          style={{
            paddingBottom: '50%',
            width: '50%',
            position: 'relative',
            height: 0
          }}
        >
          <div
            style={{
              position: 'absolute',
              top: '0',
              left: '0',
              width: '100%',
              height: '100%'
            }}
          >
            <ResponsiveContainer>
              <PieChart>
                <Tooltip />
                <Pie
                  data={this.state.totalQueueLength}
                  dataKey="value"
                  nameKey="name"
                  cx="50%"
                  cy="50%"
                  outerRadius={'100%'}
                  fill="#8884d8"
                >
                  {this.state.totalQueueLength.map((entry, index) => (
                    <Cell key={index} fill={COLORS[index % COLORS.length]} />
                  ))}
                </Pie>
              </PieChart>
            </ResponsiveContainer>
          </div>
        </div>
      </div>
    );
  }
}

export default App;
